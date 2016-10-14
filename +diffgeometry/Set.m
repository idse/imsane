classdef Set < handle_light
    % Square subset of R^n
       
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
        
        name        % name that identifies the set
        dimension   % dimension of the set
    end

    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = Set(varargin)
            % SET Create a new set
            %
            % Set(name, dimension)
            % Set(xmlnode)
            %
            % xmlnode to restore object from xml file
            
            
            if nargin == 2 
                
                this.name = varargin{1};
                this.dimension = varargin{2};
            
            % load from xml node    
            elseif nargin == 1
                
                node = varargin{1};
                this.name = char(xpathText(node, 'name'));
                this.dimension = str2num(xpathText(node, 'dimension'));
            end
        end
        
        function hasss = hasSubset(this, set)
            % hasSubset Test whether some set is contained in this set
            %
            % hasss = hasSubset(set)
            %
            % FiniteSet overloads this so whenever this function is called
            % the set is infinite.
            % A finite set of the same dimension is always a subset,
            % an infinite set only when it is the same set.
            
            if strcmp(set.name, this.name)
                hasss = true;                
            elseif isa(set, 'FiniteSet') && set.dimension == this.dimension
                hasss = true;
            else
                hasss = false;
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
            
            objectNode = docNode.createElement('Set');
            objectNode.setAttribute('subclass', class(this));
            
            elem_node = docNode.createElement('name');
            text_node = docNode.createTextNode(this.name);
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
            
            elem_node = docNode.createElement('dimension');
            text_node = docNode.createTextNode(num2str(this.dimension));
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
        end
    end
end