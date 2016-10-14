classdef SquareTensorPatch < diffgeometry.TensorPatch & diffgeometry.SquareMatPatch
    % A tensor that is also a square matrix, i.e. something whose
    % components transform as a tensor but of which we can take determinant
    % and trace etc. E.g. metric or second fundamental form (curvature) L
       
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
        
        function this = SquareTensorPatch(varargin)
            % TENSORPATCH Create new tensor patch
            %
            % TensorPatch(domain, image, def, type, chartName, atlas)
            % TensorPatch(domain, image, def, TransformableCP)
            % ComponentPatch(XMLnode, dataDir, atlas)
            %
            % The second is for subclasses to have a contructor with the
            % same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            % The third is for loading saved maps and is called from the
            % manifold2D.load
            %
            % See also Map ComponentPatch TranformableCP

            % need to call both parent constructors for multiple
            % inheritence
            % however, if SquareMatPatch is called with XMLnode, it will
            % try to also load the tiffs, and moreover from the wrong
            % directory, so we need to call it differently
            % we can't use properties already set before calling a
            % superclass constructor, so a temporary object that reads the
            % tiffs needs to be made....
            
            tmp = diffgeometry.TensorPatch(varargin{:});
            
            % parent constructor 1
            this = this@diffgeometry.TensorPatch(tmp.domain, tmp.image, tmp.phi, tmp.type, tmp.chartName, tmp.atlas);
            
            % parent constructor 2
            this = this@diffgeometry.SquareMatPatch(tmp.domain, tmp.image, tmp.phi, tmp.type);
            
            % check that def is actually a matrix
            if ~iscell(this.phi) || ~ismatrix(this.phi) || size(this.phi, 1) ~= size(this.phi, 2)
                error('Matrix definition needs to be an NxN cell array');
            end 
        end
    end
end