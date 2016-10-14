function se = strel3D(shape, size)
    % 3D version of matlabs 2d strel function
    % Implements 'sphere' and 'cube'
    %
    % strel3d(shap,size) 
    %
    % see also strel
    
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
    
    N = size;
    
    if strcmp(shape, 'sphere')
        
        se = false([2*N+1 2*N+1 2*N+1]);
        [X,Y,Z] = meshgrid(-N:N, -N:N, -N:N);
        %se(X.^2 + Y.^2 + Z.^2 <= N^2 + 1) = 1;
        se(X.^2 + Y.^2 + Z.^2 <= N^2) = 1;

    elseif strcmp(shape, 'cube')
        
        se = true([2*N+1 2*N+1 2*N+1]);
        
    else
        
        error('strel type not recognized');
    end

end