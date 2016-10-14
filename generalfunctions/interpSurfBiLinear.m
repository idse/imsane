function result = interpSurfBiLinear(inputGrids, stack, ssGrids, ss)
    % interpolate 2D data set over gridded 2D inputs
    %
    % interpSurfBiLinear(inputgrids, stack)
    %
    % inputGrids is cell array of 2D grids {X,Y}
    % stack is a 2D numeric array
    
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
    
    
    % the ssGrids are not uniformly spaced, but uniform*ss; Rescale space,
    % interpolate and then done! 
    
    %[ssXu,ssYu] = meshgrid(1:size(ssGrids{1},1),1:size(ssGrids{2},2));


    %disp('blabla');
    
    % rescale inputGrid by ss;
    inputGridsScaled = {single(inputGrids{1}/ss)+1,single(inputGrids{2}/ss)+1};


    % now perform similar to intersurftrilinear; 

    inputSize = size(inputGridsScaled{1});

    % we cannot properly interpolate on boundary pixels so we will exclude
    % them by setting to NaN
    for i = 1:2
        inputGridsScaled{i}(inputGridsScaled{i}<1) = NaN;
        inputGridsScaled{i}(inputGridsScaled{i}>size(ssGrids{i},i)) = NaN;
    end
    %
    % determine for which grid values any of the grids are NaNs
    isNaNGrids = false(inputSize);
    for i = 1:2
        isNaNGrids( isnan(inputGridsScaled{i}(:)) ) = true;
    end
    sind = [2,1];

    for i = 1:2
        % exclude NaNs from interpolation
        inputGridsScaled{i} = inputGridsScaled{i}(~isNaNGrids);
        % vector of pixel indices at which we are reading out the stack
        floorInd{i} = floor(inputGridsScaled{i}(:)) ;
        % subpixel displacement
        floorInd{i}(floorInd{i}>size(stack,sind(i))-1) = size(stack,sind(i))-1;
        floorInd{i}(floorInd{i}<0) = 0;
        eps{i} = single(inputGridsScaled{i} -  floor(inputGridsScaled{i})); 
        % the computation of eps becomes a little more complicated ...  

    end


    c3 = {};
    for i = 0:1
        for j = 0:1
            linind = sub2ind(size(stack), floorInd{2} + j,...
                floorInd{1} + i);
            c3{i + 1,j + 1} = single(stack(linind));
        end
    end
    
    
    % second linear interpolation, on 2 faces parallel to the 12-plane
    c2 = {};
    for j = 0:1
        c2{j + 1} = (1-eps{1}).*c3{1, j + 1} + eps{1}.*c3{2, j + 1};
    end

    c0 = (1-eps{2}).*c2{1} + eps{2}.*c2{2};

    %c0 = (1-eps{3}).*c1{1} + eps{3}.*c1{2};
   
    % put the vectorized result back on the grid 
    result = zeros(inputSize);
    result(~isNaNGrids) = c0;
end