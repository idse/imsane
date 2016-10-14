classdef Manifold2D < handle_light
    % 2+1 dimensional manifold: a surface in time
    % 
    % Time dependence is dealt with by making arrays of Fields indexed by
    % a time index (most of the time called ti). If the Fields representing
    % the geometry (embedding and metric) are themselves time dependent
    % the property dynamic = true and embedding and metric are also arrays, 
    % but this is not always the case.
       
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
        
        atlas            % object containing charts and transition maps
        fields           % cell array of fields on this surface
        topologicalSpace % object containting sets and intersections
        
        dynamic          % boolean: is the geometry dynamic
        timePoints       % times at which want to define time dependent fields
    end
    
    %---------------------------------------------------------------------
    % dependent properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected, Dependent = true)
        
        g               % short hand for metric tensor
        embedding       % object holding embedding of SOI in the data stack
        nTimePoints;    % number of time points
    end
    
    %---------------------------------------------------------------------
    % protected methods
    %---------------------------------------------------------------------
        
    methods (Access= protected)
        
        function fieldIndex = getFieldIndex(this, fieldName)
            % GETFIELDINDEX Get field index, return 0 if nonexistent
            %
            % fieldIndex = getFieldIndex(fieldName)
            
            for fieldIndex=1:length(this.fields)
                if strcmp(this.fields{fieldIndex}(1).name, fieldName)
                    return;
                end
            end
            fieldIndex = 0;
        end
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = Manifold2D(varargin)
            % MANIFOLD2D Create a 2D manifold (surface)
            % 
            % Manifold2D(embeddingSpace, timePoints, dynamic)
            % Manifold2D(datadir)
            %
            % embeddingSpace:   FiniteSet
            % timePoints:       vector of integer timepoints for which fields
            %                   should be defined
            % dynamic:          if false, the geometry is time independent
            % datadir:          directory holding a saved manifold
            
            % if loading the surface
            if nargin == 1
                
                dataDir = varargin{1};
                this.load(dataDir);
                
            % else creating new surface
            else
                embeddingSpace = varargin{1};
                assert(isa(embeddingSpace,'diffgeometry.Set'), 'embeddingSpace should be a set');
                
                this.timePoints = varargin{2};
                assert(isnumeric(this.timePoints), 'timePoints should be vector');
                
                this.dynamic = varargin{3};
                assert(this.dynamic==0 | this.dynamic==1, 'dynamic should be number');

                if this.dynamic
                    arrayLen = this.nTimePoints;
                else
                    arrayLen = 1;
                end

                % make topological space and atlas
                topSpaceArray(arrayLen) = diffgeometry.TopologicalSpace;

                % the atlas for each time needs to be fed a pointer to the
                % topological space at the same time
                atlasArray(arrayLen) = diffgeometry.Atlas;
                for i = 1:arrayLen
                    atlasArray(i) = diffgeometry.Atlas(topSpaceArray(i));
                end

                this.topologicalSpace = topSpaceArray;
                this.atlas = atlasArray;

                % make empty list of fields
                this.fields = {};

                % create the embedding field
                this.createField('embedding', 'diffgeometry.CoordinateMap',...
                                                    embeddingSpace, this.dynamic);

                % create the metric field
                targetSpace = diffgeometry.Set('metric', 4);
                this.createField('metric', 'diffgeometry.SquareTensorPatch',...
                                                    targetSpace, this.dynamic);
            end
        end

        %------------------------------------------------------
        % load from disc
        %------------------------------------------------------
        
        function load(this, dataDir)
            % LOAD Load saved 2D manifold
            % 
            % load(dataDir)

            % reading of actual image files happens in Map constructor
            
            DOMnode = xmlread(fullfile(dataDir, 'SOI.xml'));

            this.dynamic = logical(str2double(xpathText(DOMnode,'SOI/dynamic')));
            this.timePoints = str2num(xpathText(DOMnode, 'SOI/timePoints'));
            node = xpathNode(DOMnode, 'SOI/fields/Field[name=''embedding'']/targetSpace/Set');
            embeddingSpace = diffgeometry.FiniteSet(node);
            
            %----------------------------------------------
            % create empty embedding, topology and atlas
            %----------------------------------------------
            
            if this.dynamic
                arrayLen = this.nTimePoints;
            else
                arrayLen = 1;
            end
            
            % make topological space and atlas
            topSpaceArray(arrayLen) = diffgeometry.TopologicalSpace;
            
            % the atlas for each time needs to be fed a pointer to the
            % topological space at the same time
            atlasArray(arrayLen) = diffgeometry.Atlas;
            for i = 1:arrayLen
                atlasArray(i) = diffgeometry.Atlas(topSpaceArray(i));
            end
            
            this.topologicalSpace = topSpaceArray;
            this.atlas = atlasArray;
            
            %------------------------------------------
            debugMsg(1, 'loading topology\n');
            %------------------------------------------
            
            % static case
            if ~this.dynamic
                
                node = xpathNode(DOMnode, 'SOI/topologicalSpace/TopologicalSpace');
                this.topologicalSpace = diffgeometry.TopologicalSpace(node);

            % dynamic case
            else

                topSpaceArrayNode = xpathNode(DOMnode, 'SOI/topologicalSpace');
                topSpaceNode = topSpaceArrayNode.getFirstChild;
                
                while ~isempty(topSpaceNode)
                    
                    if strcmp(topSpaceNode.getNodeName, 'TopologicalSpace')
                        
                        t = str2num(topSpaceNode.getAttribute('time'));
                        tidx = this.tIdx(t);
                        this.topologicalSpace(tidx) = diffgeometry.TopologicalSpace(topSpaceNode);
                    end 
                    topSpaceNode = topSpaceNode.getNextSibling;
                end
            end
            
            %------------------------------------------
            debugMsg(1, 'loading atlas\n');
            %------------------------------------------
            
            % static case
            if ~this.dynamic

                node = xpathNode(DOMnode, 'SOI/atlas/Atlas');
                this.atlas = diffgeometry.Atlas(this.topologicalSpace, node, dataDir);
                
            % dynamic case
            else
                
                atlasArrayNode = xpathNode(DOMnode, 'SOI/atlas');
                atlasNode = atlasArrayNode.getFirstChild;
                
                while ~isempty(atlasNode)
                    
                    if strcmp(atlasNode.getNodeName, 'Atlas')
                        
                        t = str2num(atlasNode.getAttribute('time'));
                        tidx = this.tIdx(t);
                        this.atlas(tidx) = diffgeometry.Atlas(this.topologicalSpace, atlasNode, dataDir);
                    end 
                    atlasNode= atlasNode.getNextSibling;
                end
            end

            %----------------------------------------------
            % create empty embedding and metric
            %----------------------------------------------
            
            % make empty list of fields
            this.fields = {};

            % create the embedding field
            this.createField('embedding', 'diffgeometry.CoordinateMap',...
                                                embeddingSpace, this.dynamic);

            % create the metric field
            targetSpace = diffgeometry.Set('metric', 4);
            this.createField('metric', 'diffgeometry.SquareTensorPatch',...
                                                targetSpace, this.dynamic);
                                            
            %------------------------------------------
            debugMsg(1, 'loading fields\n');
            %------------------------------------------
            
            fieldsNode = xpathNode(DOMnode, 'SOI/fields');

            fieldNode = fieldsNode.getFirstChild;

            while ~isempty(fieldNode)

                if strcmp(fieldNode.getNodeName, 'Field')

                    t = str2num(fieldNode.getAttribute('time'));
                    fieldName = char(xpathText(fieldNode, 'name'));
                    patchClass = char(xpathText(fieldNode, 'patchClass'));

                    debugMsg(2, [fieldName ', t = ' num2str(t) '\n']);
                    
                    node = xpathNode(fieldNode, 'targetSpace/Set');
                    tspaceClass= str2func(char(node.getAttribute('subclass')));
                    targetSpace = tspaceClass(node);

                    field = this.getField(fieldName);

                    if field == 0
                        % for now, assume only embedding and metric can be static 
                        % they are constructed earlier so never here
                        fDynamic = 1; 
                        this.createField(fieldName, patchClass, targetSpace, fDynamic);
                        field = this.getField(fieldName);
                    end

                    % the atlas that needs to be passed to patch
                    % constructor
                    if this.dynamic
                        gti = this.tIdx(t);
                    else
                        gti = 1;
                    end
                    patlas = this.atlas(gti);
                    
                    % patches: contain the actual field data
                    patchesNode = xpathNode(fieldNode, 'patches');
                    patchNode = patchesNode.getFirstChild;
                    
                    while ~isempty(patchNode)

                        if strcmp(patchNode.getNodeName, 'Map')

                            patchConstructor = str2func(patchClass);
                            domName = xpathText(patchNode, 'domain/Set/name');
                            pdir = fullfile(dataDir, 'fields', fieldName, domName);
                            patch = patchConstructor(patchNode, pdir, patlas);

                            % if subsampled, make fullsize (see Atlas)
                            if patch.subSampling > 1
                                patch = patchConstructor(patch.domain, patch.image, patch.apply, patch);
                            end
                            
                            % if the data cannot be found, phi is empty
                            if ~isempty(patch.phi)
                                field(this.tIdx(t)).addPatch(patch);
                            end
                        end
                        
                        patchNode = patchNode.getNextSibling;
                    end

                    % ALSO RENAME XMLNODE FUNCTION IN SET CLASS TO SAVE
                end
                fieldNode = fieldNode.getNextSibling;
            end
        end
          
        %------------------------------------------------------
        % save
        %------------------------------------------------------
        
        function save(this, options)
            % SAVE: saves SOI metadata in SOI.xml and data in tiffs
            %
            % save(options)
            %
            % options:              struct with the following fields
            %   - dir               directory to save to
            %   - imwriteOptions    options to pass to imwrite
            
            % directory to hold the SOI
            assert(isfield(options, 'dir'), 'options.dir needs to be specified');
            dir = options.dir;
            if ~exist(dir, 'dir'), mkdir(dir); debugMsg(1, ['created dir ' dir]); end
            
            % by default, don't save metric
            if ~isfield(options, 'saveMetric'); options.saveMetric = false; end
            
            % default imwrite options use tif and deflate
            if ~isfield(options, 'imwriteOptions')
                options.imwriteOptions = {'tif', 'Compression', 'deflate'};
            end
            
            %-----------------
            % create XML
            %-----------------
            
            docNode = com.mathworks.xml.XMLUtils.createDocument('SOI');
            timeCreated = int16(datevec(now));
            docNode.getDocumentElement.setAttribute('timeStamp', num2str(timeCreated, '%d'));
            
            % dynamic
            elem_node = docNode.createElement('dynamic');
            text_node = docNode.createTextNode(num2str(this.dynamic));
            elem_node.appendChild(text_node);
            docNode.getDocumentElement.appendChild(elem_node);

            % timePoints
            elem_node = docNode.createElement('timePoints');
            text_node = docNode.createTextNode(num2str(this.timePoints));
            elem_node.appendChild(text_node);
            docNode.getDocumentElement.appendChild(elem_node);

            %-------------------------------------------
            debugMsg(1, 'Saving atlas\n');
            %-------------------------------------------
            
            atlasopts = options;
            atlasopts.dir = fullfile(dir, 'atlas'); 
            
            elem_node = docNode.createElement('atlas');
            
            for gti = 1:length(this.atlas)
                
                % save atlas
                fnamePostfix = ['T' num2str(this.timePoints(gti), '%04d')];
                
                atlasNode = this.atlas(gti).save(docNode, atlasopts, fnamePostfix);
                if this.dynamic
                    atlasNode.setAttribute('time', num2str(this.timePoints(gti)));
                end
                elem_node.appendChild(atlasNode);
            end
            
            docNode.getDocumentElement.appendChild(elem_node);

            %-------------------------------------------
            debugMsg(1, 'Saving topology\n');
            %-------------------------------------------
            
            elem_node = docNode.createElement('topologicalSpace');
            
            for gti = 1:length(this.atlas)

                % topologicalSpace
                topSpaceNode = this.topologicalSpace(gti).save(docNode);
                if this.dynamic
                    topSpaceNode.setAttribute('time', num2str(this.timePoints(gti)));
                end
                elem_node.appendChild(topSpaceNode);
            end
            
            docNode.getDocumentElement.appendChild(elem_node);
            
            %---------------------------------
            debugMsg(1, 'Saving fields\n');
            %---------------------------------
            
            % data directory structure:
            % /fieldname/patchdomain/rep/comp_i_j_k_time_t.tif
            
            % fields can have different metadata at different times because patches may
            % have different index type from raising or lowering for example
            % we are redundantly storing lots of metadata like targetspace but this
            % simply reflects the way we deal with timedependence

            fields_node = docNode.createElement('fields');
            docNode.getDocumentElement.appendChild(fields_node);
            
            % optionally exclude metric from saving
            if options.saveMetric
                fieldIndices = 1:numel(this.fields);
            else
                fieldIndices = [1 3:numel(this.fields)];
            end
            
            % for each field
            for f = fieldIndices
                
                field = this.fields{f};
                fieldName = field(1).name;
                debugMsg(2, [fieldName '\n']);
                
                % create directory
                fdir = fullfile(dir, 'fields', field(1).name);
                if ~exist(fdir, 'dir'), mkdir(fdir); end 
                
                % options for field save
                fopts = options;
                fopts.dir = fdir;
                
                % we subsample the embedding
                if strcmp(fieldName, 'embedding')
                    fopts.subsample = 2;
                end
                
                % for each time point of the field
                tmax = length(field);
                for ti = 1:tmax
                    
                    fnamePostfix = ['T' num2str(this.timePoints(ti), '%04d')];
                    fieldNode = field(ti).save(docNode, fopts, fnamePostfix);
                    
                    fieldNode.setAttribute('time', num2str(this.timePoints(ti)));
                    fields_node.appendChild(fieldNode);
                end
            end

            xmlwrite(fullfile(dir, 'SOI.xml'), docNode);
        end
       
        %------------------------------------------------------
        % fields related
        %------------------------------------------------------
        
        function field = getField(this, fieldName, varargin)
            % GETFIELD Get field, return 0 if nonexistent
            %
            % field = getField(fieldName)
            %
            % field will generally be a Field array indexed by time
            
            fieldIndex = this.getFieldIndex(fieldName);
            
            if fieldIndex>0
                field = this.fields{fieldIndex};
            else
                field = 0;
            end
        end
        
        function names = getFieldNames(this)
            % GETFIELDNAMES Return names of all fields defined
            
            nFields = length(this.fields);
            names = cell([nFields 1]);
            for i=1:nFields
                names{i} = this.fields{i}.name;
            end
        end
        
        function removeField(this, fieldName)
            % Remove a field
            %
            % removeField(fieldName)
            
            fieldIndex = this.getFieldIndex(fieldName);
            this.fields(fieldIndex) = [];
        end
        
        function createField(this, fieldName, patchType, targetSpace, dynamic)
            % CREATEFIELD Create a field (function) on the surface
            %
            % createField(fieldName, patchType, targetSpace, dynamic)
            %
            % example: the metric
            % targetSpace = diffgeometry.Set('metric', 4);
            % this.createField('metric', 'diffgeometry.TensorPatch',...
            %                                   targetSpace, this.dynamic);
            
            if dynamic, arrayLen = this.nTimePoints;
            else        arrayLen = 1;                   
            end

            % creation of a field is like creation of the atlas, every time
            % needs to be linked to the top space at the same time

            fieldArray(arrayLen) = diffgeometry.Field;
            
            for i = 1:arrayLen
                
                % if the geometry is time dependent (as opposed to the
                % field being created) the topological space is different
                % for each time
                if this.dynamic 
                    topSpace = this.topologicalSpace(i);
                else
                    topSpace = this.topologicalSpace(1);
                end
                    
                fieldArray(i) = diffgeometry.Field(fieldName, patchType,...
                                                    targetSpace, topSpace);
            end 
            
            fieldIndex = this.getFieldIndex(fieldArray(1).name);
            
            if fieldIndex == 0
                ntensors = length(this.fields);
                this.fields{ntensors + 1} = fieldArray;
            else
                disp(['field ' fieldName ' already exists, overwriting']);
                this.fields{fieldIndex} = fieldArray;
            end
        end
             
        %------------------------------------------------------
        % calculate induced metric and add to fields
        %------------------------------------------------------

        function NCalcInducedMetric(this, varargin)
            % NCALCINDUCEDMETRIC Numerical computation of the induced metric
            %
            % NCalcInducedMetric()
            % NCalcInducedMetric(t)
            %
            % For dynamic SOI, time argument needs to be provided.
            
            if this.dynamic == false
                gti = 1;
                if length(varargin) == 1
                    debugMsg(2, 'NCalcInducedMetric: static SOI: ignoring time argument'); 
                end
            else
                if length(varargin) == 1
                    gti = this.tIdx(varargin{1});
                else error('for dynamic SOI, time argument needs to be provided'); 
                end
            end

            targetSpace = this.g(1).targetSpace;
            
            % for each chart, calculate the induced metric as the Jacobian
            % of the embedding over the surface coordinates contracted with
            % itself
            for i = 1:length(this.atlas(gti).charts)

                chart = this.atlas(gti).charts{i};
                chartDom = chart.domain.name;
                J = this.embedding(gti).getPatch(chartDom).getJacobian(chart);

                gP = J.contract(1, 1, J);
                
                if this.g(gti).getPatch(gP.domain.name) == 0
                    debugMsg(2, ['NCalcInducedMetric: adding new patch ' chart.image.name '\n']);
                    
                    % 'cast' from ComponentPatch to TensorPatch
                    gP = diffgeometry.SquareTensorPatch(gP.domain, targetSpace, gP.apply,...
                                            gP.type, chart.image.name, this.atlas(gti));
                    this.g(gti).addPatch(gP);
                else
                    debugMsg(2, 'NCalcInducedMetric: patch exists -> transform');
                    patch = this.g(gti).getPatch(gP.domain.name);
                    chartInv = chart.getInverse();
                    
                    % compose with chart inverse to get metric in chart
                    gPinChart = gP.compose(chartInv);
                    patch.setTransform(chart.image.name, gPinChart.apply, gP.type);
                end
            end
        end
        
        %------------------------------------------------------
        % calculate affine connection and add to fields
        % NOT THOROUGHLY TESTED
        % CALCULATES THE RIGHT THINGS?
        % WORKS FOR DYNAMICAL GEOMETRY? 
        % CORRECT FOR STEPSIZE IN DERIVATIVES?
        %------------------------------------------------------
        
        function NCalcChristoffel(this, chartName, varargin)
            % NCALCCHRISTOFFEL Calculate the affine connection
            %
            % NCalcChristoffel(chartName)
            % NCalcChristoffel(chartName, t)
            %
            % Creates or overwrites a ChristoffelPatch Field called connection
            %
            % for dynamic SOI, time argument needs to be provided
            % chartName = [] will do all charts
            
            if this.dynamic == false
                gti = 1;
                if length(varargin) == 1
                    debugMsg(2, 'NCalcChristoffel: static SOI: ignoring time argument'); 
                end
            else
                if length(varargin) == 1
                    gti = this.tIdx(varargin{1});
                else error('for dynamic SOI, time argument needs to be provided'); 
                end
            end
            
            % create the connection field if it doesnt exist
            if this.getField('connection') == 0
                targetSpace = diffgeometry.Set('connection', 8);
                this.createField('connection', 'diffgeometry.ChristoffelPatch',...
                                                    targetSpace, this.dynamic);
            end
            Christoffel = this.getField('connection');
            
            % the charts to compute the connection for
            if ~isempty(chartName)
                ci = this.atlas.getChartIndex(chartName);
                if ci == 0
                    error('chart is not in atlas');
                end
            else
                ci = 1:length(this.atlas(gti).charts);
            end
            
            % for each chart,
            for i = ci

                chart = this.atlas(gti).charts{i};
                gP = this.g.getPatch(chart.domain.name);
                
                % derivatives of the metric
                [dgxxdx, dgxxdy] = GaussGradient(gP.cmp({1,1}), 1);
                [dgyydx, dgyydy] = GaussGradient(gP.cmp({2,2}), 1);
                [dgxydx, dgxydy] = GaussGradient(gP.cmp({1,2}), 1);
                
                % first create patches with first index lower
                def = {};
                def{1,1,1} = dgxxdx/2;
                def{1,1,2} = dgxydx - dgxxdy/2;
                def{2,2,1} = dgxydy - dgyydx/2;
                def{2,2,2} = dgyydy/2;
                def{1,2,1} = -dgxxdy/2;
                def{1,2,2} = dgyydx/2 - dgxydy;
                def{2,1,1} = def{1,2,1};
                def{2,1,2} = def{1,2,2};
                
                ChristP = diffgeometry.ChristoffelPatch(chart.domain,...
                                Christoffel.targetSpace, def, [-1 -1 -1],...
                                chart.image.name, this.atlas(gti));
                
                % now raise the last index
                ChristP.raise(3, this.g);
                            
                if Christoffel.getPatch(chart.domain.name) == 0
                    
                    debugMsg(2, ['NCalcChristoffel: adding new patch ' chart.domain.name ' \n']);
                    this.getField('connection').addPatch(ChristP);
                    
                else
                    debugMsg(2, ['NCalcChristoffel: patch ' chart.domain.name 'exists -> transform ' chart.image.name]);
                    patch = Christoffel.getPatch(chart.domain.name);
                    chartInv = chart.getInverse();
                    
                    % compose with chart inverse to get Christoffels in chart
                    ChrPinChart = ChristP.compose(chartInv);
                    patch.setTransform(chart.image.name,...
                                        ChrPinChart.apply, ChristP.type);
                end
            end
        end
        
        %------------------------------------------------------
        % Calculate the curvature
        % NOT THOROUGHLY TESTED
        % WORKS FOR DYNAMICAL GEOMETRY? 
        %------------------------------------------------------
        
        function NCalcCurvature(this, chartName, varargin)
            % NCALCCURVATURE Calculate the curvature
            %
            % NCalcCurvature(chartName)
            % NCalcCurvature(chartName, t)
            %
            % Creates or overwrites a SquareTensorPatch Field called curvature
            % The matrix L_{ij} = X^I_{ij} n^I, i.e. the derivatives
            % of the embedding with respect to the surface coordinates,
            % dotted in three dimensions into the normal. 
            % We raise one index and store L^i_j so that 
            % the Gaussian curvature is K = det(L^i_j)
            % and the mean curvature is H = tr(L^i_j)/2 
            %
            % for dynamic SOI, time argument needs to be provided
            % chartName = [] will do all charts

            if this.dynamic == false
                gti = 1;
                if length(varargin) == 1
                    debugMsg(2, 'NCalcCurvature: static SOI: ignoring time argument'); 
                end
            else
                if length(varargin) == 1
                    gti = this.tIdx(varargin{1});
                else error('for dynamic SOI, time argument needs to be provided'); 
                end
            end
            
            % create the curvature field if it doesnt exist
            if this.getField('curvature') == 0
                targetSpace = diffgeometry.Set('curvature', 4);
                this.createField('curvature', 'diffgeometry.SquareTensorPatch',...
                                                    targetSpace, this.dynamic);
            end
            curvature = this.getField('curvature');

            % the charts to compute the connection for
            if ~isempty(chartName)
                ci = this.atlas(gti).getChartIndex(chartName);
                if ci == 0
                    error('chart is not in atlas');
                end
            else
                ci = 1:length(this.atlas(gti).charts);
            end
            
            % for each chart,
            for i = ci

                % get the embedding
                chart = this.atlas(gti).charts{i};
                
                % THIS WILL NOT GIVE THE RIGHT EMBEDDING FOR
                % CYLINDER_PROPER
                X = this.embedding(gti).getPatch(chart.domain.name).apply;

                % we may want to use sigma in the future, then you want to
                % probably keep things consistent and use a metric that was
                % calculated with the same sigma, so for now just do
                % everything with sigma = 0, i.e. pixel difference gradient
                sigma = 0;
                
                % derivative of the embedding
                [Xu,Xv] = GaussGradient(X{1}, sigma);
                [Yu,Yv] = GaussGradient(X{2}, sigma);
                [Zu,Zv] = GaussGradient(X{3}, sigma);

                % second derivatives of embedding
                [Xuu,Xuv] = GaussGradient(Xu, sigma);
                [Yuu,Yuv] = GaussGradient(Yu, sigma);
                [Zuu,Zuv] = GaussGradient(Zu, sigma);
                [Xuv,Xvv] = GaussGradient(Xv, sigma);
                [Yuv,Yvv] = GaussGradient(Yv, sigma);
                [Zuv,Zvv] = GaussGradient(Zv, sigma);

                % Reshape 2D arrays into vectors and calculate normal
                Xu      =   [Xu(:) Yu(:) Zu(:)];
                Xv      =   [Xv(:) Yv(:) Zv(:)];
                normal  =   cross(Xu,Xv,2);
                normal  =   normal./repmat(sqrt(dot(normal,normal,2)), [1 3]);

                % second fundamental form : curvautre matrix
                Xuu = [Xuu(:) Yuu(:) Zuu(:)];
                Xuv = [Xuv(:) Yuv(:) Zuv(:)];
                Xvv = [Xvv(:) Yvv(:) Zvv(:)];
                Luu = dot(Xuu,normal,2);
                Luv = dot(Xuv,normal,2);
                Lvv = dot(Xvv,normal,2);

                [s,t] = size(X{1});
                Luu = reshape(Luu,s,t);
                Luv = reshape(Luv,s,t);
                Lvv = reshape(Lvv,s,t);
                
                % now correct derivative for chart stepsize
                du = chart.image.stepSize(1);
                dv = chart.image.stepSize(2);
                def = {Luu/du^2, Luv/(du*dv); Luv/(du*dv) Lvv/dv^2};
                
                curvPatch = diffgeometry.SquareTensorPatch(chart.domain,...
                                curvature.targetSpace, def, [-1 -1],...
                                chart.image.name, this.atlas);

                % raise the index of the curvature tensor so H and K are
                % simply trace and determinant
                curvPatch.raise(1, this.g(gti));
                            
                if curvature(gti).getPatch(chart.domain.name) == 0
                    
                    debugMsg(2, ['NCalcCurvature: adding new patch ' chart.domain.name ' \n']);
                    curvature(gti).addPatch(curvPatch);
                    
                else
                    debugMsg(2, ['NCalcCurvature: patch ' chart.domain.name ' exists -> transform ' chart.image.name '\n']);
                    patch = curvature(gti).getPatch(chart.domain.name);
                    
                    % if the chart is not the defining one for the patch
                    % compose with chart inverse to get curvature in chart
                    if ~strcmp(chart.image.name, patch.chartName)
                        chartInv = chart.getInverse();
                        curvPatchinChart = curvPatch.compose(chartInv);
                    else
                        curvPatchinChart = patch;
                    end
                    
                    patch.setTransform(chart.image.name,...
                                        curvPatchinChart.apply, curvPatchinChart.type);
                end
            end
        end

        %------------------------------------------------------
        % convert time to timePoint array index
        %------------------------------------------------------
        
        function idx = tIdx(this, t)
            % TIDX Convert time to timePoint array index
            %
            % idx = tIdx(t)
            % 
            % e.g. if timePoints = [0 50], tIdx(50) returns 2
            
            for idx = 1:length(this.timePoints)
                if this.timePoints(idx) == t, return; end
            end
        end
                
        %------------------------------------------------------
        % getters for dependent properties
        %------------------------------------------------------
        
        function g = get.g(this)
            % return the metric tensor
            g = this.getField('metric');
        end
        
        function embedding = get.embedding(this)
            % return the metric tensor
            embedding = this.getField('embedding');
        end
        
        function nTimePoints = get.nTimePoints(this)
            nTimePoints = numel(this.timePoints);
        end
	end
end
