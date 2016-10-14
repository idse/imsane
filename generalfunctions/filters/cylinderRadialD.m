function dstackdr = cylinderRadialD(stack, sigma, order, varargin)
    % Gaussian derivative with respect to cylindrical radial coordinate
    %
    % cylinderRadialD(stack, sigma, order)
    % cylinderRadialD(stack, sigma, order, zind)
    % cylinderRadialD(stack, sigma, order, zind, center)
    %
    % zind is the dimension representing the z-axis if the data is 3D, default z=3
    % center is a pair x,y defining the cylinder axis

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

    D = ndims(stack);
    
    % direction of cylinder axis
    if nargin > 4, zind = varargin{1};
    else zind =3; end
        
    % stack size
    s = size(stack);
    
    % for the 2D case
    if numel(s) == 2
        s(3) = 1;
    end
    
    % transverse position of cylinder axis (c=center)
    if nargin == 5
        c = [varargin{2}(2), varargin{2}(1), 0];
    else
        c = round(s/2);
    end

    % third element of idxPerm is z dimension, first two are the remaining
    % dimensions
    idxPerm = [setxor(1:D, zind), zind];
    
    dstackdr = stack;

    for i = 1:order
        
        dyx = {dstackdr, dstackdr};

        % take xy derivatives of stack
        for j = 1:2
            for k = 1:2
                % first derivative if j == k
                %dyx{i} = convn(dyx{i}, GaussD(sigma, i==j, j), 'same');
                dyx{j} = myconvn(dyx{j}, GaussD(sigma, j==k, idxPerm(k)));
            end
        end

        % make a meshgrid for the radial coordinate        
        X = {};
        [X{1}, X{2}, X{3}] = ndgrid(-c(1):(s(1)-c(1))-1,...
                            -c(2):(s(2)-c(2))-1, -c(3):(s(3)-c(3))-1);
        
        Y = X{idxPerm(1)};
        X = X{idxPerm(2)};
        R = sqrt(X.^2 + Y.^2);

        % take the nth radial derivative
        % partial_r = (x/r)partial_x + (y/r)partial_y

        dstackdr = (Y./R).*dyx{1} + (X./R).*dyx{2};
    end
end