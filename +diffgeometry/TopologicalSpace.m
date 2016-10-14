classdef TopologicalSpace < handle_light
    % Collection of sets on which maps such as chart and embedding are
    % defined.
    % Roughly corresponds to the concept of a topological space but really
    % just a bookkeeping device to keep track of domains in the surface on
    % which maps can be defined.
       
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
    
    properties (SetAccess = private)
        
        sets;       % cell array of set objects
        nsets;      % number of sets
        intersects; % matrix of which sets intersect each other
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = TopologicalSpace(varargin)
            % TOPOLOGICALSPACE Create a topologicalSpace
            %
            % TopologicalSpace()
            % TopologicalSpace(xmlnode)
            %
            % xmlnode to restore object from xml file
            
            if nargin == 0
                
                this.nsets = 0;
                this.intersects = [];
                this.sets = {};
                
            elseif nargin == 1
                
                node = varargin{1};
                
                this.nsets = str2num(xpathText(node, 'nsets'));
                this.intersects = str2num(xpathText(node, 'intersects'));
                
                this.sets = {};
                for i = 1:this.nsets
                    setnode = xpathNode(node, ['sets/Set[' num2str(i) ']']);
                    this.sets{i} = diffgeometry.FiniteSet(setnode);
                end
            end
        end
        
        % ------------------------------------------------------
        % add / remove sets
        % ------------------------------------------------------
        
        function ID = addSet(this, set, intersections) 
            % ADDSET Add a set to the topology
            %
            % ID = addSet(this, set, intersections) 
            %
            % set:              FiniteSet
            % intersections:    vector of other sets that set intersects
            % ID:               index of set in TopologicalSpace.sets
            
            
            % check that it is not in there already
            if this.getSetID(set.name) ~= 0
                disp(['set ' set.name ' is already in the topological space'...
                                                     ', overwriting']);
                this.sets{this.getSetID(set.name)} = set;
                return;
            end
            
            nsets = this.nsets;
            ID = nsets + 1;
            
            % add set to sets array
            this.sets{ID} = set;
            
            % obtain IDs of intersecting sets
            isectIDs = zeros([1 length(intersections)], 'uint8');
            for i = 1: length(intersections)
                isectIDs(i) = this.getSetID(intersections{i});
                if isectIDs(i) == 0
                    error(['intersecting set ' intersections{i}...
                        ' not part of topological space!'])
                end
            end
            
            % update intersection matrix
            tmp = false(nsets + 1);
            tmp(1:nsets, 1:nsets) = this.intersects;
            tmp(ID, isectIDs) = true;
            tmp(isectIDs, ID) = true;
            % the trivial intersection is always there
            tmp(ID, ID) = true;
            
            this.intersects = tmp;
            
            % update this.nsets
            this.nsets = nsets + 1;
        end
        
        function setID = getSetID(this, setName) 
            % GETSETID Return setID (index in .sets) or zero if nonexistent.
            %
            % setID = getSetID(setName)
            
            if ~ischar(setName)
                error('argument should be a setName string');
            end
            
            for setID=1:this.nsets
               if strcmp(this.sets{setID}.name, setName)
                   return;
               end
            end
            setID = 0;
        end
        
        function intersections = getIntersections(this, setName) 
            % GETINTERSECTIONS Return names of sets that have a non-empty 
            % intersection with a particular set
            %
            % intersections = getIntersections(setName)
            %
            % intersections:    cell array of names
            
            ID = this.getSetID(setName);
            intersections = {};
            if ID > 0
                for i = find(this.intersects(ID, :))
                    intersections = [intersections, this.sets{i}.name];
                end
            else
                error('set not in topology');
            end
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
            
            objectNode = docNode.createElement('TopologicalSpace');
            
            elem_node = docNode.createElement('nsets');
            text_node = docNode.createTextNode(num2str(this.nsets));
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
            
            elem_node = docNode.createElement('intersects');
            text_node = docNode.createTextNode(vec2str(uint8(this.intersects)));
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
            
            elem_node = docNode.createElement('sets');
            for i = 1:length(this.sets)
                subelem_node = this.sets{i}.save(docNode);
                elem_node.appendChild(subelem_node);
            end
            objectNode.appendChild(elem_node);
        end
    end
end