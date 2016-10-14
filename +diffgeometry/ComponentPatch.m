classdef ComponentPatch < diffgeometry.Map
    % Multi-component object that can be contracted as a tensor 
    % (and therefore raised and lowered) 
    % but includes non-tensors like Jacobian and connection.

    % TODO: save time + memory by storing and using index permutation symmetries
       
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
        
        % type - vector of index types, -1 lower, 0 neutral, 1 upper
        % e.g. e.g. metric: [-1 -1], Jacobian [-1 0], Christoffel [1 -1 -1]
        type;
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------

    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = ComponentPatch(varargin)
            % COMPONENTPATCH Create a component patch
            %
            % ComponentPatch(domain, image, def, type)
            % ComponentPatch(domain, image, def, CP)
            % ComponentPatch(XMLnode, dataDir, atlas)
            %
            % The second is for subclasses to have a contructor with the
            % same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            % The third is for loading saved maps and is called from the
            % manifold2D.load
            %
            % See end of code of ComponentPatch.contractOther for example
            % of calling the constructor with an object.
            %
            % See also Map
            
            if nargin == 4
                
                if isa(varargin{4}, 'diffgeometry.ComponentPatch')
                    CP = varargin{4};   % ComponentPatch from which to take type
                    type = CP.type;
                else
                    type = varargin{4};
                end
                
                varargin = varargin(1:3);
                
            elseif nargin == 3
                
                XMLnode = varargin{1};
                type = str2num(xpathText(XMLnode, 'type'));

            else
                error('wrong number of arguments');
            end

            % call parent constructor
            this = this@diffgeometry.Map(varargin{:});
            
            if isempty(this.phi), return; end
            
            % a vector should be a column
            % test should be after calling map to avoid dealing with the
            % subsampled structure input
            if numel(size(this.phi)) == 2 && all(size(this.phi) == [1 2])
                debugMsg(1, 'ComponentPatch: transposing definition to get  column vector\n');
                this.phi = this.phi';
            end
            
            % scalar, map should have one element
            % vector, map should be a column of length 2
            % something else should at least have a type specified for each
            % dimension
            % remember, type can be 0, i.e. non-tensor index
            if  isempty(type) && numel(this.phi)~=1 ||... 
                numel(type) == 1 && type(1) ~= 0 && ~all(size(this.phi)==[2 1]) ||...
            	numel(type) > 1 && numel(type) ~= numel(size(this.phi))
                error([ 'type length ' num2str(numel(type))...
                        ' does not match CP definition dimension '...
                        num2str(numel(size(this.phi)))]);
            else
                this.type = type;
            end
        end
        
        % ------------------------------------------------------
        % get/set components
        % ------------------------------------------------------
        
        function comp = cmp(this, index, varargin)
            % CMP Return component over input if image matches required
            % chartname.
            % If no input given, return component over entire domain.
            %
            % cmp(index)
            % cmp(index, input)
            %
            % this function calls applyComponent
            % See also applyComponent
            
            if nargin == 3
                input = varargin{1};
            else
                input = this.domain;
            end
            
            comp = this.applyComponent(index, input);
        end
        
        % ------------------------------------------------------
        % contract indices
        % ------------------------------------------------------
        
        function newCP = contract(this, index, otherIndex, varargin)
            % CONTRACT Contract with other ComponentPatch or self
            % 
            % contract(this, index, otherIndex)         with self
            % contract(this, index, otherIndex, other)  with other
            %
            % See also contractSelf, contractOther

            if ~isnumeric(index) || ~isnumeric(otherIndex)
                error('first two arguments of contract should be indices');
            end
            
            if nargin == 4
                % contraction with other tensor
                otherCP = varargin{1};
                newCP = this.contractOther(index, otherIndex, otherCP);
            else
                debugMsg(1, 'Contract with self\n');
                newCP = this.contractSelf(index, otherIndex);
            end    
        end
                
        % ------------------------------------------------------
        % raise and lower indices
        % ------------------------------------------------------
        
        function lower(this, index, g)
            % LOWER Lower an index with metric g
            %
            % lower(index, g)
            %
            % metric has to be passed because patch doesn't know about time
            % dependence so we need to feed it the metric at the right time
            
            if index > length(this.type)
                error('lowering index is greater than number of indices');
            end
            if this.type(index) ~= 1
                disp('index is already lower or cannot be lowered, doing nothing');
                return;
            end
            
            % if the component patch is defined on an index set this should
            % return nonzero
            if g.getPatchIndex(this.domain.name) > 0
                
                gPatch = g.getPatch(this.domain.name);
            
            % if it is defined on a chart, we have to find which index set
            % corresponds to that chart to get the right patch
            % because passing the atlas complicates things, we do it
            % without
            else
                
                idx = [];
                for i = 1:length(g.patches)
                    
                    idx = find(strcmp(g.patches{i}.availableReps, this.domain.name));
                    
                    if ~isempty(idx)
                        chartName = g.patches{i}.availableReps{idx};
                        gPatch = g.patches{i}.getTransform(chartName);
                    end
                end
            end
    
            if ~exist('gPatch', 'var')
                error('It appears that the metric is not defined on the domain of this patch.');
            end
            
            lowered = this.contract(index, 1, gPatch);
            
            % contract append indices of the second tensor at the end, for
            % lowering we want the lowered index to be in the same place so
            % we have to rearrange things
            % unless we're talking about a vector
            lphi = lowered.phi;
            ltype = lowered.type;
            
            if length(ltype) > 1
                % rearranged order of indices
                endidx = length(ltype);
                idxRe = [1:index-1 endidx index:endidx-1];

                % rearrange definition dimensions
                phiRe = permute(lphi, idxRe);
                typeRe = ltype(idxRe);

                % overwrite current definition
                this.type = typeRe;
                this.phi = phiRe;
            else
                this.type = ltype;
                this.phi = lphi;
            end
        end
        
        function raise(this, index, g)
            % RAISE Raise the index with the metric g
            %
            % raise(index, g)
            %
            % metric has to be passed because patch doesn't know about time
            % dependence so we need to feed it the metric at the right time
            
            if index > length(this.type)
                error('raising index is greater than number of indices');
            end
            if this.type(index) ~= -1
                disp('index is already upper or cannot be raised, doing nothing');
                return;
            end
            
            % if the component patch is defined on an index set this should
            % return nonzero
            if g.getPatchIndex(this.domain.name) > 0
                
                gP = g.getPatch(this.domain.name);
            
            % if it is defined on a chart, we have to find which index set
            % corresponds to that chart to get the right patch
            % because passing the atlas complicates things, we do it
            % without
            else
                
                idx = [];
                for i = 1:length(g.patches)
                    
                    idx = find(strcmp(g.patches{i}.availableReps, this.domain.name));
                    
                    if ~isempty(idx)
                        chartName = g.patches{i}.availableReps{idx};
                        gP = g.patches{i}.getTransform(chartName);
                    end
                end
            end
            
            if ~exist('gP', 'var')
                error('It appears that the metric is not defined on the domain of this patch.');
            end

            % we compute the inverse metric every time
            gPInv = gP.inverse();
            raised = this.contract(index, 1, gPInv);
            
            % contract append indices of the second tensor at the end, for
            % raising we want the raised index to be in the same place so
            % we have to rearrange things
            % unless it is a vector
            rphi = raised.phi;
            rtype = raised.type;
            
            if length(rtype) > 1
                % rearranged order of indices
                endidx = length(rtype);
                idxRe = [1:index-1 endidx index:endidx-1];

                % rearrange definition dimensions
                phiRe = permute(rphi, idxRe);
                typeRe = rtype(idxRe);

                % overwrite current definition
                this.type = typeRe;
                this.phi = phiRe;
            else
                % overwrite current definition
                this.type = rtype;
                this.phi = rphi;
            end
        end
        
        % ------------------------------------------------------
        % change image name
        % ------------------------------------------------------
        
        function setImageName(this, name)
            % SETIMAGENAME Set name of ComponentPatch image
            %
            % setImageName(name)
            
            if isa(this, 'diffgeometry.FiniteSet');
                error('TODO: implement setImageName for FiniteSet image');
            else
                this.image = diffgeometry.Set(name, this.image.dimension);
            end
        end
        
        %------------------------------------------------------
        % XML
        %------------------------------------------------------
        
        function objectNode = save(this, docNode, options, varargin)
            % SAVE: saves TransformableCP data to tiff and returns metadata
            % as XML node
            % overloads Map.save to add index type to XML
            %
            % objectNode = save(docNode, options)
            % objectNode = save(docNode, options, filenamePostfix)
            %
            % docNode:              document node to which objectNode should belong
            % objectNode:           node representing the object
            % options:              struct with the following fields
            %   - dir               directory to save to
            %   - imwriteOptions    options to pass to imwrite
            %
            % tiff filename is of form 'cmp_i_j_filenamePostfix.tif'
            % filenamePostfix can be used to append time for example
            %
            % See also Map.save
            
            objectNode = save@diffgeometry.Map(this, docNode, options, varargin{:});

            % type
            elemNode = docNode.createElement('type');
            text_node = docNode.createTextNode(num2str(this.type));
            elemNode.appendChild(text_node);
            objectNode.appendChild(elemNode);
        end
    end
    
    %---------------------------------------------------------------------
    % protected methods
    %---------------------------------------------------------------------

    methods (Access = protected)
        
        % ------------------------------------------------------
        % contract with self
        % ------------------------------------------------------
        
        function newCP = contractSelf(this, index, otherIndex)
            % CONTRACTSELF Contract two ComponentPatch indices
            % 
            % newCP = contractSelf(index, otherIndex)
            %
            % newCP object is of same class as this
            
            % order the indices for convenience
            inds = sort([index otherIndex], 'ascend');
            index = inds(1);
            otherIndex = inds(2);
            
            % see if we are contracting contravariant with covariant
            if this.type(index) ~= -this.type(otherIndex)
                error('can only contract indices of opposite type');
            end
      
            % the new index type vector and the index ranges
            notIndex = setdiff(1:length(this.phiSize), [index otherIndex]);

            newType = this.type(notIndex);
            newIndRanges  = this.phiSize(notIndex);
               
            % the number of indices of the resulting CP
            nIndNew = length(newType);
            
            % initialize the new component cell array
            if nIndNew > 1
                newGrids = cell(newIndRanges);

            % but the syntax is cell([d 1]) for a vector
            elseif nIndNew == 1
                newGrids = cell([newIndRanges 1]);

            % and for a scalar it is just one component
            elseif nIndNew == 0
                newGrids = cell(1);
                newType = []; % this is expected for a scalar
            end
            
            for i = 1:numel(newGrids)
                newGrids{i} = zeros([this.domain.gridSize(2), this.domain.gridSize(1)]);
            end
            
            % all the index values of the newCP to loop over
            indexVals = this.indexValues(newIndRanges);
            
            % for each combination of index values of the new CP
            for p = 1: size(indexVals,1)
                
                newInd = num2cell(indexVals(p,:));
                
                % sum the product of this and other tensor over the
                % contracted index
                % e.g. if we contract this_{kili}
                % then newCP_{kl} = sum_i this_{kili}

                for i = 1:this.phiSize(index)
                    thisInd = num2cell([indexVals(p, 1:index-1) i ...
                                        indexVals(p, index:otherIndex-2) i
                                        indexVals(p, otherIndex:nIndNew)]);
                         
                    newGrids{newInd{:}} = newGrids{newInd{:}} + this.cmp(thisInd);
                end
            end
            
            % TODO: nameing could be improved
            imageName = [this.image.name '_' num2str(index) '-' num2str(otherIndex)];
            image = diffgeometry.Set(imageName, numel(newGrids));
            
            % construct newCP object, of same class as this if constructor
            % arguments match ComponentPatch (otherwise override this)
            CPConstructor = str2func(class(this)); 
            newCP = CPConstructor(this.domain, image, newGrids, newType);
        end
        
        
        % ------------------------------------------------------
        % contract with other
        % ------------------------------------------------------
        
        function newCP = contractOther(this, index, otherIndex, otherCP)
            % CONTRACTOTHER Contract with other component patch
            % 
            % newCP = contractOther(index, otherIndex, otherCP)
            %
            % newCP object is of same class as this
            
            % see if we are contracting contravariant with covariant
            if this.type(index) ~= -otherCP.type(otherIndex)
                error('can only contract indices of opposite type');
            end
            % see that we are contracting components in the same
            % coordinates
            if ~strcmp(this.domain.name, otherCP.domain.name)
               error(['Contracting component patch on"' this.domain.name...
                   '" with other CP on "' otherCP.domain.name '" is impossible']);
            end
            
            % generate an instructive message
            debugMsg(1, ['Contracting [' num2str(this.type(index)) '] index '...
                num2str(index) ' of ' this.image.name ' with ['...
                num2str(otherCP.type(otherIndex))...
                '] index ' num2str(otherIndex') ' of ' otherCP.image.name '\n']);

            % the new index type vector and the index ranges
            notIndex = setdiff(1:length(this.type), index);
            notOtherIndex = setdiff(1:length(otherCP.type), otherIndex);

            newType = [this.type(notIndex) otherCP.type(notOtherIndex)];
            newIndRanges  = [this.phiSize(notIndex) otherCP.phiSize(notOtherIndex)];
            
            % the number of indices of this CP and the resulting CP
            nIndThis = length(this.type);
            nIndNew = numel(newIndRanges);
            
            % this is just because of the stupid syntax of functions like
            % cell which take cell(3) to mean a 3x3, instead 3x1
            if nIndNew > 1
                newGrids = cell(newIndRanges);

            % but the syntax is cell([d 1]) for a vector
            elseif nIndNew == 1
                newGrids = cell([newIndRanges 1]);

            % and for a scalar it is just one component
            elseif nIndNew == 0
                newGrids = cell(1);
                newType = []; % this is expected for a scalar
            end
            
            % for subsampled grids, the contracted grid should be the
            % smaller of two fields being contracted
            if this.subSampling > 1 || otherCP.subSampling > 1
                
                domGrids = this.domain.makeGrids;
                
                if this.subSampling > otherCP.subSampling
                    
                    newGridSize = size(this.phi{1});
                    
                    ssZero = this.ssZero;
                    ss = this.subSampling;
                    
                    x = domGrids{1}(ssZero:ss:end, ssZero:ss:end);
                    y = domGrids{2}(ssZero:ss:end, ssZero:ss:end);
                    
                    thisGrids = this.phi;
                    otherGrids = otherCP.apply({x,y});
                    
                    newSS = this.subSampling;
                    newSS0 = this.ssZero;
                else
                    newGridSize = size(otherCP.phi{1});
                    
                    ssZero = otherCP.ssZero;
                    ss = otherCP.subSampling;
                    
                    x = domGrids{1}(ssZero:ss:end, ssZero:ss:end);
                    y = domGrids{2}(ssZero:ss:end, ssZero:ss:end);
                    
                    thisGrids = this.apply({x,y});
                    otherGrids = otherCP.phi;
                    
                    newSS = otherCP.subSampling;
                    newSS0 = otherCP.ssZero;
                end
                
            % no subsampling is easy
            else
                newGridSize = [this.domain.gridSize(2), this.domain.gridSize(1)];
                thisGrids = this.apply;
                otherGrids = otherCP.apply;
                newSS = 1;
                newSS0 = 1;
            end
            
            % initialize the new component cell array
            for i = 1:numel(newGrids)
                newGrids{i} = zeros(newGridSize);
            end

            % all the index values of the newCP to loop over
            indexVals = this.indexValues(newIndRanges);
            
            % for each combination of index values of the new CP
            for p = 1: size(indexVals,1)
                
                newInd = num2cell(indexVals(p,:));
                
                % sum the product of this and other tensor over the
                % contracted index
                % e.g. if we contract this_{kil} other^i_{mn}
                % then newCP_{klmn} = sum_i this_{kil}*other^i_{mn}
                
                for i = 1:this.phiSize(index)
                    
                    % e.g. thisInd = kil, so [newInd(1) i newInd(2)] 
                    thisInd = num2cell([indexVals(p, 1:index-1) i ...
                                           indexVals(p, index:nIndThis-1)]);
                    % e.g. otherInd = imn, so [i newInd(3:4)] 
                    otherInd = num2cell([indexVals(p, nIndThis:nIndThis + otherIndex-2)... 
                             i indexVals(p, nIndThis+otherIndex-1:nIndNew)]);
                      
                    newGrids{newInd{:}} = newGrids{newInd{:}} + ...
                                    thisGrids{thisInd{:}}.*otherGrids{otherInd{:}};
                end
            end
            
            % TODO: naming could be improved
            imageName = [this.image.name '_' otherCP.image.name '_'...
                                    num2str(index) '-' num2str(otherIndex)];
            image = diffgeometry.Set(imageName, numel(newGrids));
            
            % construct newCP object
            % yeah, this is not pretty: if I contract a tensor I want to keep
            % a tensor (or TransformableCP), but this has a different
            % constructor, because of Map composition there is a universal
            % constructor call as well where I pass the object itself to
            % define the remaining properties, but this object has the
            % wrong index type so I work around that.
            % It works
            CPConstructor = str2func(class(this)); 
            tmpType = this.type;
            
            % definition is struct for subsampling
            if newSS > 1
                newDef = struct('grids', {newGrids}, ...
                    'subSampling', newSS, 'ssZero', newSS0);
            else
                newDef = newGrids;
            end
            
            this.type = newType;
            newCP = CPConstructor(this.domain, image, newDef, this);
            this.type = tmpType;
        end
        
    end
end