classdef ChristoffelPatch < diffgeometry.TransformableCP
    % TransformableCP that transforms like the affine connection
       
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
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = ChristoffelPatch(domain, image, def, varargin)
            % CHRISTOFFELPATCH Create new patch for the connection
            %
            % ChristoffelPatch(domain, image, def, type, chartName, atlas)
            % ChristoffelPatch(domain, image, def, TransformableCP)
            %
            % The latter is for subclasses of Map to have a contructor with 
            % the same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            %
            % See also Map ComponentPatch TranformableCP
            
            % parent constructor
            this = this@diffgeometry.TransformableCP(domain, image, def, varargin{:});
        end
    end
    
    methods 
        
        function newCP = calculateTransform(this, chartName)
            % CALCULATETRANSFORM Not implemented yet
            disp('Christoffel transformation not implemented yet');
        end
    end
end