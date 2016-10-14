classdef FiniteSet < diffgeometry.Set
    % Square subset of R^n that can be turned into a grid
       
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
        
        % boundary - cell array containing set boundaries, 
        % For example, if set in R^2 is (a,b)x(c,d) where (a,b) are intervals 
        % on the real line then boundary is {[a b], [c d]}
        boundary        
        stepSize	% vector with step sizes
        gridSize    % gridSize 
        
        % the information in gridSize and boundary is redundant, but
        % rounding errors in boundary can lead to problems, so we only
        % use the lower boundary, stepSize and gridSize to reload,
        % guaranteeing gridSize is always the same
    end

    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = FiniteSet(varargin)
            % FINITESET Create finite set
            %
            % FiniteSet(name, boundary, stepSize)
            % FiniteSet(xmlnode)
            %
            % xmlnode to restore object from xml file
           
            % for some bullshit reason matlab doesn't allow conditional
            % superclass constructor calls, hence the args construction
            % stepSize length = dimension
            if nargin==3, superArg = [varargin(1) length(varargin{3})];
            else superArg = varargin(1); end
            
            this = this@diffgeometry.Set(superArg{:});

            % properties from arguments
            if nargin==3
            
                this.name = varargin{1};
                this.boundary = varargin{2};
                this.stepSize = varargin{3};

                assert(length(this.boundary) == length(this.stepSize),...
                    ['boundary dimension ' num2str(length(this.boundary)) ' and stepSize length ' num2str(length(this.stepSize)) ' do not match']);
            
                % gridSize
                this.gridSize = zeros([1 this.dimension]);
                for i = 1:this.dimension
                    this.gridSize(i) = floor((this.boundary{i}(2) - this.boundary{i}(1))/this.stepSize(i) + 1);
                end
                
            % properties from xml
            elseif nargin == 1
                
                node = varargin{1};
                this.stepSize = str2num(xpathText(node, 'stepSize'));
                
                this.boundary = {};
                for i = 1:this.dimension
                    path = ['boundary/direction[@index=' num2str(i) ']'];
                    this.boundary{i} = str2num(xpathText(node, path));
                end
                
                % load gridSize
                this.gridSize = str2num(xpathText(node, 'gridSize'));
                
                % now enforce consistency of boundary with gridSize, because
                % this can be lost by rounding errors
                for i = 1:length(this.boundary)
                    this.boundary{i}(2) = this.boundary{i}(1) + (this.gridSize(i) - 1)*this.stepSize(i);
                end
            end
            
            % check some boundary properties
            assert(~isempty(this.boundary), 'argument boundary is empty');
            
            for i = 1:length(this.boundary)
                assert(this.boundary{i}(1) <= this.boundary{i}(2),...
                    ['boundary{' num2str(i) '} has lower limit > upper limit']);
            end
        end

        
        %------------------------------------------------------
        % grids, stepsize
        %------------------------------------------------------
        
        function grids = makeGrids(this)
            % MAKEGRIDS Return meshgrid of set values
            %
            % makeGrids()

            linear = {};
            
            for i=1:this.dimension    
                range = this.boundary{i};
                du = this.stepSize(i);
                linear{i} = single(range(1):du:range(2));
            end
            
            if this.dimension == 2

                [u1, u2] = meshgrid(linear{1}, linear{2});
                grids = {u1, u2};
                
            elseif this.dimension == 3
                
                [u1, u2, u3] = meshgrid(linear{1}, linear{2}, linear{3});
                grids = {u1, u2, u3};
                
            else
                error('dimension should be 2 or 3');
            end
        end
        
        function handles = makeHandles(this)
            % MAKEHANDLES Return handle for (linear) function from indices to set value
            %
            % handles = makeHandles()
            %
            % handles:  cell array of handles indexed by dimension
            %
            % This gives an analytic expression relating indices in the Set
            % grid to values in the set grid.
            
            handles = cell([1 this.dimension]);
            for i=1: this.dimension
                handles{i} = @(x) this.stepSize(i)*(x{i}-1) + this.boundary{i}(1);
            end
        end
        
        function handles = makeInverseHandles(this)
            % MAKEINVERSEHANDLES Handle for (linear) function from set value to indices 
            % 
            % handles = makeInverseHandles()
            %
            % handles:  cell array of handles indexed by dimension
            %
            % This gives an analytic expression relating values in the set 
            % grid to indices in the set grid.
            % 
            % See also FiniteSet.makeHandles
            
            handles = cell([1 this.dimension]);
            for i=1: this.dimension
                handles{i} = @(x) (x{i} - this.boundary{i}(1))/this.stepSize(i) + 1;
            end
        end
        
        %------------------------------------------------------
        % subset
        %------------------------------------------------------
        
        function hasss = hasSubset(this, set)
            % HASSUBSET Test whether some set is contained in this set
            %
            % hasss = hasSubset(set)
            
            % an infinite set or set of different dimension is no subset
            if ~isa(set, 'diffgeometry.FiniteSet') || this.dimension ~= set.dimension
                hasss = false;  
                return;

            % otherwise see if the set is fully contained in this one
            else
                hasss = true;
                for i = 1:this.dimension
                    hasss = hasss & (set.boundary{i}(1) >= this.boundary{i}(1));
                    hasss = hasss & (set.boundary{i}(2) <= this.boundary{i}(2));
                end
            end
        end
        
        function indices = getSubsetIndices(this, set)
            % GETSUBSETINDICES Return linear indices for a subset of the set
            %
            % indices = getSubsetIndices(set)
            
            if ~ this.hasSubset(set)
                error('indices can only be returned for a subset');
            end
            
            indices = cell(size(this.boundary));
            
            for i = 1:this.dimension
                a = floor((set.boundary{i}(1) - this.boundary{i}(1))/this.stepSize(i) + 1);
                b = floor((set.boundary{i}(2) - this.boundary{i}(1))/this.stepSize(i) + 1);
                indices{i} = a:b;
            end
        end
        
                
        %------------------------------------------------------
        % index grid
        %------------------------------------------------------
     
        function indexSet = makeIndexSet(this)
            % MAKEINDEXSET Create a set of grid indices of this set
            %
            % indexSet = makeIndexSet()
            
            indexName = [this.name,'_index'];
            indexBoundary  = {[1 this.gridSize(1)], [1 this.gridSize(2)]};
            indexStepSize  = [1 1];
            indexSet  = diffgeometry.FiniteSet(indexName, indexBoundary, indexStepSize);
        end
        
        %------------------------------------------------------
        % XML
        %------------------------------------------------------
     
        function objectNode = save(this, docNode)
            % XMLNODE Create an XML node representing instance of this class
            %
            % objectNode = XMLnode(docNode)
            % 
            % docNode : document node to which objectNode should belong
            % objectNode : node representing the object
            
            objectNode = save@diffgeometry.Set(this, docNode);
            
            elem_node = docNode.createElement('stepSize');
            text_node = docNode.createTextNode(vec2str(this.stepSize));
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
            
            % the information in gridSize and boundary is redundant, see
            % comment in properties
            elem_node = docNode.createElement('gridSize');
            text_node = docNode.createTextNode(vec2str(this.gridSize));
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
            
            elem_node = docNode.createElement('boundary');
            
            for i = 1:this.dimension
                
                subelem_node = docNode.createElement('direction');
                subelem_node.setAttribute('index', num2str(i));
                text_node = docNode.createTextNode(vec2str(this.boundary{i}));
                subelem_node.appendChild(text_node);
                elem_node.appendChild(subelem_node);
            end

            objectNode.appendChild(elem_node);
        end
    end
end