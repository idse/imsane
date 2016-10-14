function [dx, dy] = GaussGradient(im, sigma)
    % GAUSSGRADIENT Like gradient but using Gaussian derivative
    %
    % GaussGradient(im, sigma)
    %
    % for sigma == 0 this just calls gradient
    
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
    
    % make sure 1D case is along x
    if size(im,2) == 1
        im = im';
    end
    
    if sigma > 0
        
        % 1D case
        if size(im,1) == 1
            
            dx = myconvn(im, GaussD(sigma, 1, 2));
            nanmask = isnan(dx);
            dx(nanmask) = 0;
            dy = [];
        
        % 2D case
        else
            dx = myconvn(im, GaussD(sigma, 1, 2));
            dy = myconvn(im, GaussD(sigma, 1, 1));

            nanmask = isnan(dx) | isnan(dy);
            dx(nanmask) = 0;
            dy(nanmask) = 0;
        end
    else
        
        % 1D case
        if size(im,1) == 1
            dx = gradient(im);
            dy = [];
        % 2D case
        else
            [dx, dy] = gradient(im);
        end
    end
end