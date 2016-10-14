classdef Map < handle_light
    % Take elements in domain set to elements in image set
       
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2015 Idse Heemskerk and Sebastian Streichan
    %
    % This file is part of ImSAnE.
    % 
    % ImSAnE is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % ImSAnE is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with ImSAnE.  If not, see <http://www.gnu.org/licenses/>.
    %
    % We make an explicit exception giving permission to link this code
    % with GPL-incompatible matlab facilities.
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected)
        
        domain;             % set on which the map is defined
        image;              % image set of map
        subSampling = 1;    % subsampling factor: imss = im(ssZero:subSampling:end);  
        ssZero = 1;         % subsampling offset
        
        phi;            % definition: cell array of grids or function handles
        % phi not to be confused with an azimuthal angle, phi is a fairly
        % standard name for a map, e.g. a chart is a pair (U, phi) where U
        % is a subset (here U = domain) and phi maps domain -> image
    end
    
    properties (Dependent = true)
        
        representation; % type of internal representation: analytic or numeric
        phiSize;        % size of the definition cell array
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = Map(varargin)
            % Map Create a map between two sets
            %
            % Map(domain, image, definition)
            % Map(domain, image, definition, map)
            % Map(XMLnode, dataDir, atlas)
            %
            % The second is for subclasses to have a contructor with the
            % same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            % The third is for loading saved maps and is called from the
            % manifold2D.load
            %
            % Definition can be:
            % - a cell array of grids for numeric rep
            % - a cell array of function handles for analytic rep
            % - a structure with properties:
            %   - grids: cell array of subsampled grids
            %   - subSampling: subsampling factor
            %   - ssZero: subsampling zero
            %   such that grids{i} = fullGrids{i}(ssZero:subSampling:end);

            if isa(varargin{1}, 'diffgeometry.Set')
               
                this.domain = varargin{1};
                this.image = varargin{2};
                definition = varargin{3};
                
                % if definition is a struct
                if isa(definition, 'struct')

                    if numel(definition) > 1
                        error('subsampled struct cannot be array, maybe pass grids as {grids}?');
                    end

                    if isfield(definition, 'subSampling')

                        this.subSampling = definition.subSampling;  
                        this.ssZero = definition.ssZero;  

                        debugMsg(2, ['Map: subSampling '...
                            num2str(this.subSampling) ' from definition struct\n']);
                    else
                        error('subsampled structure is missing field subSampling');
                    end 

                    if isfield(definition, 'grids')
                        definition = definition.grids;
                    else
                        error('if definition is struct it should contain a field grids');
                    end
                end
            
            elseif isa(varargin{3}, 'diffgeometry.Atlas')
                
                node = varargin{1};
                dataDir = varargin{2};

                domNode = xpathNode(node, 'domain/Set');
                domClass = str2func(char(domNode.getAttribute('subclass')));
                this.domain = domClass(domNode);
                
                imNode = xpathNode(node, 'image/Set');
                imClass= str2func(char(imNode.getAttribute('subclass')));
                this.image = imClass(imNode);
                
                fnamePostfix = xpathText(node, 'filenamePostfix');
                
                % read map definition from tiffs
                phiSize = str2num(xpathText(node, 'phiSize'));
                grids = cell(phiSize);
                
                indexVals = this.indexValues(phiSize);
                for i = 1:size(indexVals, 1)
                    
                    fname = ['cmp' num2str(indexVals(i,:), '_%d') fnamePostfix];
                    
                    % three different file extensions we support
                    % or maybe this should be stored?
                    if exist(fullfile(dataDir, [fname '.tif']),'file')
                        im = imread(fullfile(dataDir, [fname '.tif']));
                    
                    elseif exist(fullfile(dataDir, [fname '.jp2']),'file')
                        im = imread(fullfile(dataDir, [fname '.jp2']));
                        
                    elseif exist(fullfile(dataDir, [fname '.jpg']),'file')
                        im = imread(fullfile(dataDir, [fname '.jpg']));
                        
                    elseif exist(fullfile(dataDir, [fname '.png']),'file')
                        im = imread(fullfile(dataDir, [fname '.png']));
                    else
                        im = [];
                        warning(['Map: file not found: ' fname ' creating empty object']);
                        return;
                    end
                    
                    idxCell = num2cell(indexVals(i,:));
                    grids{idxCell{:}} = im;
                end
                
                definition = grids;
                
                ssNode = xpathNode(node, 'subSampling');
                if ~isempty(ssNode)
                    this.subSampling = str2num(ssNode.getTextContent);
                    this.ssZero = str2num(xpathText(node, 'ssZero'));
                end
            else
                error('arguments of wrong type');
            end

            % check map definition
            assert(iscell(definition), 'map is specified by a cell array of grids or function handles');
            
            assert( numel(definition) == this.image.dimension,...
                ['number of grid cell array elements ' num2str(numel(definition))...
                 ' does not match image dimension ' num2str(this.image.dimension)]);
            
            
            if isa(definition{1}, 'function_handle')
                
                assert(isa(this.domain, 'diffgeometry.Set'), 'domain should be a Set object (not domain.name)');
                this.phi = definition;
            
            % if numeric, check size and domain type
            elseif isnumeric(definition{1}) || islogical(definition{1})
                
                if ~isa(this.domain, 'diffgeometry.FiniteSet')
                   error('numeric map definition requires domain to be FiniteSet');
                end
                
                % check that definition size matches domain size up to
                % subsampling
                for i = 1:numel(definition)
                    
                    if this.subSampling == 1
                        domainSize = this.domain.gridSize;
                    else
                        % see subsample()
                        domainSize = floor((this.domain.gridSize -this.ssZero + 1)/this.subSampling) + 1;
                    end
                    
                    domainSize(1:2) = domainSize([2 1]); % matlab sucks yx -> xy
                    
                    if any(domainSize ~= size(definition{i})) 
                        error(['numerical map definition size ' num2str(size(definition{i}))...
                            ' should match domain ' this.domain.name...
                            ' size (yxz, ceil((size-ss0+1)/subsampling)) ' num2str(domainSize)]);
                    end
                end
                
                this.phi = definition;
                
            else
                error('definition cell array should contain grids or function handles');
            end            
        end

        %------------------------------------------------------
        % apply map to set
        %------------------------------------------------------
        
        function result = apply(this, varargin)
            % APPLY Let map act on finite set or list of input values
            %
            % result = APPLY(S)
            % result = APPLY({p1,..,pn}) 
            % result = APPLY()
            %
            % S:        finite set in the domain 
            % {p1,..}:  cell array containing arrays of input values for
            %           each domain dimension
            % result:   cell array of output values
            % If no input is given, apply to entire domain.
            %
            % This really just calls applyComponent in a loop.
            % Therefore, when a single component is needed, it is more
            % efficient to call applyComponent.'
            %
            % See also applyComponent

            % deal with the subsampled case
            if this.subSampling ~= 1
                fullMap = this.upsample();
                result = fullMap.apply;
                return
            end
            
            % initialize result cell array
            result = cell(this.phiSize);
            
            % fill result cell array
            for i=1:numel(result)
                result{i} = this.applyComponent({i}, varargin{:});
            end
        end
        
        function result = applyComponent(this, ind, varargin)
            % applyComponent Act with a component of a map on finite set or 
            % input values
            %
            % result = applyComponent(ind, S)
            % result = applyComponent(ind, {p1,..,pn})
            % result = applyComponent(ind)
            % 
            % ind :     index values specified as cell array, e.g. {1,2}
            % S:        finite set in the domain 
            % {p1,..}:  cell array containing arrays of input values for
            %           each domain dimension
            % result:   grid of output values
            % If no second argument is given, apply to entire domain.
            %
            % Two kinds of input are useful because while arrays of
            % input values are most general, it requires generating them.
            % On the other hand a set only specifies boundaries and
            % stepsize: if an input set just calls part of the map grid 
            % we dont have to grid the input or compute anything
            
            if ~iscell(ind) || numel(ind{1})~=1
                error('indices should be specified as cell array {i,j,..}');
            end
            
            if nargin == 3
                input = varargin{1};
            else
                input = this.domain;
            end
            
            % use the right definition to generate result:
            
            debugMsg(3, ['Map.applyComponent ' num2str(cell2mat(ind)) ': '...
                this.domain.name ' -> ' this.image.name ';  ' this.representation ';']);
            
            tID = tic;
            
            %-------------------------------
            % if the input is a finite set
            %-------------------------------
            if isa(input, 'diffgeometry.FiniteSet')
                
                if ~this.domain.hasSubset(input)
                    error('Map can only act on a set contained in its domain');
                end
                
                % if the input is a set and stepsize matches the domain stepsize
                % then we can just return part of the grid, if the grid exists
                if all(input.stepSize == this.domain.stepSize)...
                                && strcmp(this.representation, 'numeric')...
                                && this.subSampling == 1;
                    
                    % if called for the entire domain, just return the entire
                    % grid
                    if strcmp(input.name, this.domain.name)
                        result = this.phi{ind{:}};                    

                    % else call part of the grid
                    else
                        range = this.domain.getSubsetIndices(input);

                        % matlab sucks! (first index, second coordinate)

                        if this.domain.dimension==2
                            result = this.phi{ind{:}}(range{2}, range{1}); 
                        elseif this.domain.dimension==3
                            result = this.phi{ind{:}}(range{2}, range{1}, range{3});
                        else
                            error('dimension is not 2 or 3?')
                        end
                    end
                % if we can't just return part of the grid, we have to
                % compute the values of the input set and then do whatever
                % we would do with a set of input grids
                % this includes analytic and subsampled maps
                else
                    inputGrids = input.makeGrids();
                    result = this.applyComponent(ind, inputGrids);
                end
            
            %-------------------------------
            % if the input is a cell array
            %-------------------------------
            elseif iscell(input)

                if numel(input) ~= this.domain.dimension
                   error(['Number of input cell array elements '...
                       num2str(numel(input))...
                       ' should match map domain dimension '...
                       num2str(this.domain.dimension)]);
                end
                
                % if this map is analytic
                if strcmp(this.representation, 'analytic')

                    result = this.phi{ind{:}}( input );
                    
                % interpolation depends both on the domain dimensions
                % we use interp2 and interp3 because they use little memory
                % and are faster than triscatteredinterp if only called a
                % few times on the same map
                elseif strcmp(this.representation, 'numeric')
                  
                    if this.domain.dimension == 2
                        
                        domainGrids = this.domain.makeGrids();

                        % interp needs conversion to single or double
                        result = interp2(...
                            domainGrids{1}, domainGrids{2},...
                            single(this.phi{ind{:}}),...   
                            input{1}, input{2}, 'linear');
                        
                    elseif this.domain.dimension == 3
                        
                        result = interpSurfTriLinear(input, this.phi{ind{:}});

                    else
                        error('dimension is not 2 or 3?');
                    end
                end            
                    
            else
                error('Map has to act on a FiniteSet or cell array of value arrays');
            end

            debugMsg(3, ['dt = ' num2str(toc(tID)) ' sec\n']);
            
            % finally, cast result type back to original output type
            if isnumeric(this.phi{1})
                result = cast(result, class(this.phi{1}));
            end
        end
        
        %------------------------------------------------------
        % compose map
        %------------------------------------------------------
        
        function composition = compose(this, mp)
            % COMPOSE Act on other map and return composed map
            %
            % composition = COMPOSE(mp)
            %
            % mp is the map to be composed with
            % composition is the result
            %
            % As functions: if m:B->C, mp:A->B, where m denotes this map, 
            % then composition = m(mp):A->C.
            % The map type of composition is therefore the type of m.

            debugMsg(1, 'Map.compose(): ');
            
            if ~isa(mp, 'diffgeometry.Map')
                error('maps can only be composed with other maps');
            end
            if ~strcmp(this.domain.name, mp.image.name)
               warning(['mp.image.name "' mp.image.name...
                   '" does not match this.domain.name "' this.domain.name '"']);
            end
            
            % numeric composition
            if strcmp(mp.representation, 'numeric') ||...
                    strcmp(this.representation, 'numeric')
                
                debugMsg(2, 'numeric\n');
                
                imp = mp.apply; 
                compDef = this.apply(imp);

            % analytic composition 
            elseif strcmp(mp.representation, 'analytic') &&...
                        strcmp(this.representation, 'analytic')
                
                debugMsg(2, 'analytic\n');
                compDef = cell(this.phiSize);

                for i = 1:numel(this.phi)

                    if numel(mp.phi) == 1
                        compDef{i} = @(u) this.phi{i}({mp.phi{1}(u)});
                    elseif numel(mp.phi) == 2
                        compDef{i} = @(u) this.phi{i}({mp.phi{1}(u), mp.phi{2}(u)});
                    elseif numel(mp.phi) == 3
                        compDef{i} = @(u) this.phi{i}({mp.phi{1}(u), mp.phi{2}(u), mp.phi{3}(u)});
                    elseif numel(mp.phi) == 4
                        compDef{i} = @(u) this.phi{i}({mp.phi{1}(u), mp.phi{2}(u), mp.phi{3}(u), mp.phi{4}(u)});
                    else
                        error('many components, expand the code');
                    end
                end
            end
            
            % the composed map has the type of the outer (this) map
            compConstructor = str2func(class(this));
            composition = compConstructor(mp.domain, this.image, compDef, this);
        end
                
        %------------------------------------------------------
        % compose map with inverse of other map
        %------------------------------------------------------

        function mMpInv = composeInverse(this, mp)
            % composeInverse Compose with inverse of other map
            % 
            % mMpInv = composeInverse(mp)
            %
            % mp:       invertibe map
            % mMpInv:   composition
            % The only invertible maps we ever compose with are coordinate
            % maps so mp has to be a coordinate map.
            %
            % As functions: if m:B->C, mp:C->A, where m denotes this map, 
            % then mMpInv = mp(m):C->A.
            % The map type of composition is therefore the type of mp.
            %
            % This function exists because it is numerically much more
            % efficient than computing the inverse separately. 
            
            if ~isa(mp, 'diffgeometry.CoordinateMap')
                disp('you have to compose with a coordinate map');
            end
            
            if ~strcmp(mp.domain.name, this.domain.name)
                error(['domain of input ' mp.domain.name...
                    ' does not match domain of map ' this.domain.name]);
            end
            
            % grids for this map
            phi = this.apply(this.domain);
            
            % grid of primed map
            phiPrimeGrids = mp.apply(mp.domain);
            
            % grid of primed image space is domain grid for this map
            grids = mp.image.makeGrids(); 
                        
            % composing this with the inverse of mp means
            % that we this over the mp image space
            % TriScatteredInterp needs double inputs...
            interpolant = {};
            phiPhipInv = {};
            
            for d = 1:numel(phi)
                interpolant{d} = TriScatteredInterp(double(phiPrimeGrids{1}(:)),...
                    double(phiPrimeGrids{2}(:)), double(phi{d}(:)));
                phiPhipInv{d} = interpolant{d}(double(grids{1}), double(grids{2}));
            end
            
            mMpInv = diffgeometry.CoordinateMap(mp.image, this.image, phiPhipInv);
        end
        
        %------------------------------------------------------
        % getters of dependent vars
        %------------------------------------------------------
        
        function rep = get.representation(this)
            
            if isa(this.phi{1}, 'function_handle') 
                rep = 'analytic';
                
            elseif isnumeric(this.phi{1}) || islogical(this.phi{1})
                rep = 'numeric';
            end
        end
        
        function phiSize = get.phiSize(this)
            
            phiSize = size(this.phi);
        end
 
        % ------------------------------------------------------
        % component index values 
        % ------------------------------------------------------
        
        function indVals = indexValues(this, varargin)
            % INDEXVALUES Get all the index values for a Map recursively
            %
            % indVals = indexValues()
            % indVals = indexValues(indexRanges)
            %
            % Needed to loop over all components of a Map, for example
            % in contract in ComponentPatch to assign values to all 
            % components of the contraction.
            %
            % Example: if we have 2 indices in 2 dimensions:
            % indexRanges = [2 2]
            % indVals = [1 1; 2 1; 1 2; 2 2]
            % we can then call this.phi{indVals(r,:)} to get a component
            %
            % if indexRanges is not passed, phiSize is used for indexRanges

            if nargin == 2
                indexRanges = varargin{1};
            else
                indexRanges = this.phiSize;
            end
            
            nIndices = length(indexRanges);
            
            if nIndices == 0
                indVals = 1;
                
            elseif nIndices == 1
                indVals = (1:indexRanges(1))';
                
            else
                % recursively call this function
                values = this.indexValues(indexRanges(1:end-1));
                
                % repeat as many times as the next index takes values
                leftblock = repmat(values, indexRanges(end), 1);
                prerightcolumn = repmat((1:indexRanges(end))', 1, size(values,1))';
                rightcolumn = reshape(prerightcolumn, numel(prerightcolumn), 1);
                indVals = [leftblock, rightcolumn];
            end
        end

        %------------------------------------------------------
        % subsample and upsample
        %------------------------------------------------------
        
        function fullMap = upsample(this)
            % create fullsize map from subsampled map
            %
            % fullMap = upsample()
            
            % for details on te code, see subsample
            
            % trivial case
            if this.subSampling == 1 
                fullMap = this;
                return
            end
            
            % full size
            xSize = this.domain.gridSize(1);
            ySize = this.domain.gridSize(2);
            
            ssGrids = this.phi;
            
            % trim last row of pixels for even sized grids
            if ~rem(xSize,2)
                ssXRange = 1:size(ssGrids{1},2)-1;
            else
                ssXRange = 1:size(ssGrids{1},2);
            end
            if ~rem(ySize,2)
                ssYRange = 1:size(ssGrids{1},1)-1;
            else
                ssYRange = 1:size(ssGrids{1},1);
            end
            
            for i = 1:numel(ssGrids)
                ssGrids{i} = ssGrids{i}(:,ssXRange);
                ssGrids{i} = ssGrids{i}(ssYRange,:);
            end

            % make domain grids and subsample them for interpolation
            domainGrids = this.domain.makeGrids();
            ss0 = this.ssZero;
            ss = this.subSampling;
            ssDomainGrids{1} = domainGrids{1}(ss0:ss:end, ss0:ss:end);
            ssDomainGrids{2} = domainGrids{2}(ss0:ss:end, ss0:ss:end);
            
            % interp needs conversion to single or double
            fullGrids = {};
            for i = 1:numel(ssGrids)
                %fullGrids{i} = interp2(...
                %    ssDomainGrids{1}, ssDomainGrids{2},...
                %    single(ssGrids{i}),...   
                %    domainGrids{1}, domainGrids{2}, 'linear');
                
                % try to interpolate using the home made function
                % interpSurfBiLinear.m
                %ssGrids{1}
                %fullGrids{i} = interpSurfBiLinear(inputGrids, ssGrids{i}, ssGrids, ss);
                
                inputGrids = {domainGrids{1},domainGrids{2}};
                fullGrids{i} = fastLinearInterp(inputGrids,ssGrids{i},ss);
                %imshow(fullGrids{i},[])
                %drawnow
                %pause(.1)
            end
            
            % restore the row / column that couldn't be interpolated
            if ~rem(xSize,2)
                for i = 1:numel(fullGrids)
                    
                    edge = interp1( domainGrids{2}(ss0:ss:end, end),...
                                    this.phi{i}(ssYRange,end),...
                                    domainGrids{2}(:, end) );
                    fullGrids{i}(:,end) = edge;
                end
            end
            if ~rem(ySize,2)
                for i = 1:numel(fullGrids)
                    
                    edge = interp1( domainGrids{1}(end,ss0:ss:end),...
                                    this.phi{i}(end,ssXRange),...
                                    domainGrids{1}(end,:) );
                    fullGrids{i}(end,:) = edge;
                    fullGrids{i}(end,end) = this.phi{i}(end,end);
                end
            end
            
            constructor = str2func(class(this));
            fullMap = constructor(this.domain, this.image, fullGrids, this);
        end
        
        function ssMap = subsample(this, ssFactor)
            % create subsampled version of Map
            % 
            % ssMap = subsample(ssFactor)
            %
            % ssFactor: integer subsampling factor
            
            % basically
            
            grids = this.apply;
            
            [ySize, xSize] = size(grids{1});
            
            yIdx = 1:ssFactor:ySize;
            xIdx = 1:ssFactor:xSize;
            
            if ~rem(ySize,2); yIdx = [yIdx, ySize]; end
            if ~rem(xSize,2); xIdx = [xIdx, xSize]; end
            
            ssGrids = cell(this.phiSize);
            for i= 1: numel(ssGrids)
                ssGrids{i} = grids{i}(yIdx, xIdx);
            end
            def = struct('grids', {ssGrids}, 'subSampling', ssFactor, 'ssZero', 1);
            
            ssConstructor = str2func(class(this));
            ssMap = ssConstructor(this.domain, this.image, def, this);
        end
        
        %------------------------------------------------------
        % save
        %------------------------------------------------------
        
        function objectNode = save(this, docNode, options, varargin)
            % SAVE Saves Map data to tiff and returns metadata as XML node
            % metadata is returned by XMLnode and linked into SOI.xml
            %
            % objectNode = save(docNode, options)
            % objectNode = save(docNode, options, filenamePostfix)
            %
            % docNode:              document node to which objectNode should belong
            % objectNode:           node representing the object
            % options:              struct with the following fields
            %   - dir               directory to save to
            %   - imwriteOptions    options to pass to imwrite
            %   - make8bit          rescale dynamic range and make 8bit
            % 
            % tiff filename is of form 'cmp_i_j_filenamePostfix.tif'
            % filenamePostfix can be used to append time for example
            
            if nargin == 4
                fnamePostfix = ['_' varargin{1}];
            else
                fnamePostfix = [];
            end
           
            assert(isfield(options, 'dir'), 'options.dir needs to be specified');
            dir = options.dir;
            
            if isfield(options, 'imwriteOptions')
                imwriteOptions = options.imwriteOptions;
            else
                imwriteOptions = {'jp2', 'Mode', 'lossless'};
            end
            
            if isfield(options, 'make8bit')
                make8bit = options.make8bit;
            else
                make8bit = false;
            end
            
            %-------------------------
            % save data to image file
            %-------------------------
            
            % if dir doesn't exist, create
            if ~exist(dir, 'dir'), mkdir(dir); disp(['created dir ' dir]); end
            
            indexVals = this.indexValues;
            
            % loop over all components
            for i = 1:size(indexVals,1)
                
                idxcell = num2cell(indexVals(i,:));
                
                % for numeric map, save definition directly because
                % applyComponent interpolates subsampled grids
                comp = this.phi{idxcell{:}};
                
                if isa(comp, 'function_handle')
                    debugMsg(1,'Saving analytic map as numeric.\n');
                    comp = this.applyComponent(idxcell);
                end
                
                % save floats as tiff
                % integers as jpeg
                filename = ['cmp' num2str(indexVals(i,:), '_%d') fnamePostfix];
                filename = fullfile(dir, filename);
                
                if isfloat(comp)
                    extension = 'tif';
                    saveFloatTiff(comp, [filename  '.' extension]);
                else
                    if make8bit
                        comp = mat2gray(comp);
                        comp = uint8(255*comp);
                    end
                    formatInfo = imformats(imwriteOptions{1});
                    extension = formatInfo.ext{1};
                    imwrite(comp, [filename  '.' extension], imwriteOptions{:});
                end
            end
            
            %---------------------
            % make XML node
            %---------------------
            
            objectNode = docNode.createElement('Map');
            objectNode.setAttribute('subclass', class(this));
           
            % domain
            elem_node = docNode.createElement('domain');
            subelem_node = this.domain.save(docNode);
            elem_node.appendChild(subelem_node);
            objectNode.appendChild(elem_node);
            
            % image
            elem_node = docNode.createElement('image');
            subelem_node = this.image.save(docNode);
            elem_node.appendChild(subelem_node);
            objectNode.appendChild(elem_node);
            
            % phiSize
            elem_node = docNode.createElement('phiSize');
            text_node = docNode.createTextNode(num2str(this.phiSize));
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
            
            % only store subSampling if it is not one
            if this.subSampling > 1
                
                elemNode = docNode.createElement('subSampling');
                text_node = docNode.createTextNode(num2str(this.subSampling));
                elemNode.appendChild(text_node);
                objectNode.appendChild(elemNode);
                
                elemNode = docNode.createElement('ssZero');
                text_node = docNode.createTextNode(num2str(this.ssZero));
                elemNode.appendChild(text_node);
                objectNode.appendChild(elemNode);
            end
            
            % filenamePostfix
            elem_node = docNode.createElement('filenamePostfix');
            text_node = docNode.createTextNode(fnamePostfix);
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
        end
    end
end