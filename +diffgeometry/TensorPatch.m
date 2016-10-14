classdef TensorPatch < diffgeometry.TransformableCP
    % TransformableCP that transforms as a tensor
       
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
        
        function this = TensorPatch(varargin)
            % TENSORPATCH Create new tensor patch
            %
            % TensorPatch(domain, image, def, type, chartName, atlas)
            % TensorPatch(domain, image, def, TransformableCP)
            % ComponentPatch(XMLnode, dataDir)
            %
            % The second is for subclasses to have a contructor with the
            % same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            % The third is for loading saved maps and is called from the
            % manifold2D.load
            %
            % See also Map ComponentPatch TranformableCP
            
            % parent constructor
            this = this@diffgeometry.TransformableCP(varargin{:});
        end
    end
    
    methods (Sealed = true)
         
        
        function newCP = calculateTransform(this, chartName)
            % CALCULATETRANSFORM Not implemented yet
            
            error('Tensor transformation not implemented yet');
        end
    end
end