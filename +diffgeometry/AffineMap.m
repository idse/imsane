classdef AffineMap < diffgeometry.CoordinateMap
    % Linear CoordinateMap defined by matrix in homogeneous coordinates
   
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
        
        matrixRep;  % matrix in homogeneous coordinates
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = AffineMap(domain, image, definition, varargin)
            % AFFINEMAP Create a map between two sets of coordinates
            %
            % AffineMap(domain, image, definition)
            % AffineMap(domain, image, definition, map)
            %
            % Definition should be homogeneous matrix, D+1 by D+1.
            % If it is not CoordinateMap constructor is called.
            %
            % See also CoordinateMap.CoordinateMap
            
            if isnumeric(definition)
                
                D = domain.dimension;
                
                % if input is numeric it must be a matrix encoding an affine
                % transformation
                if  ~ismatrix(definition) ||...
                    size(definition, 1) ~= domain.dimension + 1 ||...
                    size(definition, 2) ~= domain.dimension + 1
                    
                    error(['If map definition is numeric array it must be '...
                        'a matrix encoding an affine transformation.'...
                        'Otherwise input cell array of grids or function handles']);
                
                elseif definition(D+1, 1:D) ~= 0
                    
                    error('bottom row of homogeneous matrix should be 00..01');
                end
                
                % if matrix definition seems right, convert it to functions
                
                a = definition(1:D, D+1);         % translation
                M = definition(1:D, 1:D);         % rotation etc

                functionalDef = cell([1 D]);

                if D == 2
                    for i = 1:2
                        functionalDef{i} = @(u)( (M(i,1)*u{1} + M(i,2)*u{2}...
                            + a(i)) );
                    end
                    
                elseif D == 3
                    for i = 1:3
                        functionalDef{i} = @(u)( (M(i,1)*u{1} + M(i,2)*u{2}...
                            + M(i,3)*u{3} + a(i)) );
                    end
                end
                
            else
                debugMsg(2, ['AffineMap constructor called with non-matrix argument,'...
                    'trying CoordinateMap constructor\n']);
                functionalDef = definition;
            end
            
            % parent constructor
            this = this@diffgeometry.CoordinateMap(domain, image, functionalDef, varargin{:});
            
            % set matrix rep
            this.matrixRep = definition;
        end
        
        %------------------------------------------------------
        % inverse
        %------------------------------------------------------
        
        function invMap = getInverse(this)
            % GETINVERSE Construct inverse map by inversion of matrix representation
            %
            % mInv = getInverse()
            
            invMap = diffgeometry.AffineMap(this.image, this.domain, inv(this.matrixRep));
        end
        
    end

end