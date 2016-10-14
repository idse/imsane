classdef SquareMatPatch < diffgeometry.ComponentPatch
    % ComponentPatch that is a square matrix
    % We can then calculate trace, determinant, matrix inverse etc.
       
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
        
        function this = SquareMatPatch(domain, image, def, varargin)
            % SQUAREMATPATCH Create square matrix patch
            %
            % SquareMatPatch(domain, image, def, type)
            % SquareMatPatch(domain, image, def, CP)
            %
            % The latter is for subclasses of Map to have a contructor with 
            % the same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            %
            % See also Map ComponentPatch
            
            % parent constructor
            this = this@diffgeometry.ComponentPatch(domain, image, def, varargin{:});
            
            % check that this is actually a matrix
            if ~iscell(this.phi) || ~ismatrix(this.phi) || size(this.phi, 1) ~= size(this.phi, 2)
                error('Matrix definition needs to be an NxN cell array');
            end 
        end

        % ------------------------------------------------------
        % matrix operations
        % ------------------------------------------------------
        
        function newCP = trace(this)
            % TRACE Calculate trace of a matrix patch
            %
            % newCP = trace()
            %
            % newCP : componentPatch representing trace (scalar)
            
            N = this.phiSize(1);
            
            traceGrid = this.applyComponent({1,1});
            for i = 2:N
                traceGrid = traceGrid + this.applyComponent({i,i});
            end
            
            type = []; % scalar
            name = [this.image.name '_trace'];
            image = diffgeometry.Set(name, 1);
            newCP = diffgeometry.ComponentPatch(this.domain, image,...
                                                        {traceGrid}, type);
        end
        
        function newCP = determinant(this)
            % DETERMINANT Calculate the determinant of a matrix patch
            %
            % newCP = determinant()
            %
            % newCP : componentPatch representing determinant (scalar)
            
            if this.phiSize(1)~=2 
                error('determinant has only been implemented for 2x2 matrices');
            end
            
            grids = this.apply;
            
            detGrid = grids{1}.*grids{4} - grids{2}.*grids{3};
            
            type = []; % scalar
            name = [this.image.name '_determinant'];
            image = diffgeometry.Set(name, 1);
            newCP = diffgeometry.ComponentPatch(this.domain, image,...
                                                        {detGrid}, type);
        end
        
        function newSMP = inverse(this)
            % INVERSE Inverse matrix patch
            %
            % newSMP = inverse()
            %
            % newSMP: SquareMatPatch with inverse matrix components
            %
            % NOTE: not the inverse of this map, which goes from some set
            % in the topology to a space of matrices and is not invertible
            
            if this.phiSize(1)~=2 
                error('inverse has only been implemented for 2x2 matrices');
            end
            
            det = this.determinant().apply{1};
            detAlmostZero = det < 10^(-10);
            if any(detAlmostZero)
                warning(['the determinant of matrix patch is close to zero'...
                    'somewhere, setting to NaN']);
                det(detAlmostZero) = NaN;
            end 
            
            grids = this.apply;
                
            inverseGrids = {grids{4}./det, -grids{3}./det;
                            -grids{2}./det, grids{1}./det};
            
            debugMsg(1, 'TODO: check matrix inverse index type\n');
            type = fliplr(-this.type);
            name = [this.image.name '_inverse'];
            image = diffgeometry.Set(name, 2^2);
            newSMP = diffgeometry.SquareMatPatch(this.domain, image,...
                                                        inverseGrids, type);
        end
       
        function newSMP = multiply(this, otherSMP)
            % MULTIPLY Matrix multiplication of matrix patches
            % 
            % newSMP = multiply(otherSMP)
            
            % multiply this matrix with another
            newSMP = this.contract(2, 1, otherSMP);
        end
        
        function [lminor, lmajor, v] = diagonalize(this)
            % DIAGONALIZE Calculate eigenvalues and major axis direction 
            %
            % [lminor, lmajor, v] = diagonalize()
            % 
            % v:                unit vector field for the major eigenvector
            % lminor, lmajor:   eigenvalue fields
            % returned as grids of the internal size, so not interpolated
            % to domain size for subsampled fields
            
            % not calling this.trace because for subsampled grids that
            % interpolates unnecessarily here
            tr = this.phi{1,1};
            for i = 2:this.phiSize(1);
                tr = tr + this.phi{i,i};
            end
            
            % same for det
            det = this.phi{1}.*this.phi{4} - this.phi{2}.*this.phi{3};
            
            lmajor = (tr + sqrt(tr.^2 - 4*det))/2;
            lminor = (tr - sqrt(tr.^2 - 4*det))/2;
            
            % (v1, v2) is the major eigenvector
            % alpha = -b/(a-lambda)
            % v2 = alpha v1
            % then normalize
            alpha = -(this.phi{1,1} - lmajor)./this.phi{1,2};
            normv = sqrt(ones(size(alpha)) + alpha.^2);
            v1 = 1./normv;
            v2 = alpha./normv;
            v = {v1,v2};
        end
    end
end