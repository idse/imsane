function yidx = pointMatch(X, Y)
    % find vertex closest point in Y for each point in X
    % 
    % yidx = pointMatch(X, Y)
    %
    % X:    Nx3 or Nx2 matrix of positions
    % Y:    Mx3 or Nx2 matrix of positions
    %
    % yidx: Nx1 matrix of indices into Y
    
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

    yidx = zeros([size(X,1) 1]);

    if size(X,2) == 3
        for i = 1:size(X,1)
            [~,idx] = min(sqrt((Y(:,1) - X(i,1)).^2 + (Y(:,2) - X(i,2)).^2 ...
                                + (Y(:,3) - X(i,3)).^2));
            yidx(i) = idx;
        end
        
    elseif size(X,2) == 2
        for i = 1:size(X,1)
            [~,idx] = min(sqrt((Y(:,1) - X(i,1)).^2 + (Y(:,2) - X(i,2)).^2));
            yidx(i) = idx;
        end
    end
end