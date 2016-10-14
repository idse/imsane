function result = interpSurfTriLinear(inputGrids, stack)
    % interpolate 3D data set over gridded 2D inputs
    %
    % interpSurfTriLinear(inputgrids, stack)
    %
    % inputGrids is cell array of 2D grids {X,Y,Z}
    % stack is a 3D numeric array
    
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
    
    inputSize = size(inputGrids{1});
    
    % we cannot properly interpolate on boundary pixels so we will exclude
    % them by setting to NaN
    for i = 1:3
        inputGrids{i}(1,:) = NaN;
        inputGrids{i}(end,:) = NaN;
        inputGrids{i}(:,1) = NaN;
        inputGrids{i}(:,end) = NaN;
    end
    
    % determine for which grid values any of the grids are NaNs
    isNaNGrids = false(inputSize);
    for i = 1:3
        isNaNGrids( isnan(inputGrids{i}(:)) ) = true;
    end
    sind = [2,1,3];
    for i = 1:3
        % exclude NaNs from interpolation
        inputGrids{i} = inputGrids{i}(~isNaNGrids);
        % vector of pixel indices at which we are reading out the stack
        floorInd{i} = uint32( floor(inputGrids{i}(:)) );
        % subpixel displacement
        floorInd{i}(floorInd{i}>size(stack,sind(i))-1) = size(stack,sind(i))-1;
        floorInd{i}(floorInd{i}<2) = 2;
        eps{i} = single(inputGrids{i} -  floor(inputGrids{i}));  
        
    end

    % first linear interpolation, on the 4 edges along direction 1
    c3 = {};
    for i = 0:1
        for j = 0:1
            for k = 0:1
                linind = sub2ind(size(stack), floorInd{2} + j,...
                    floorInd{1} + i, floorInd{3} + k);
                c3{i + 1,j + 1,k +1} = single(stack(linind));
            end
        end
    end
    
    % second linear interpolation, on 2 faces parallel to the 12-plane
    c2 = {};
    for j = 0:1
        for k = 0:1
            c2{j + 1, k + 1} = (1-eps{1}).*c3{1, j + 1,k +1} + eps{1}.*c3{2, j + 1,k +1};
        end
    end

    % third linear interpolation, in the volume
    c1 = {};
    for k = 0:1
        c1{k + 1} = (1-eps{2}).*c2{1, k +1} + eps{2}.*c2{2, k +1};
    end

    c0 = (1-eps{3}).*c1{1} + eps{3}.*c1{2};
   
    % put the vectorized result back on the grid 
    result = zeros(inputSize);
    result(~isNaNGrids) = c0;
end