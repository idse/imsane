function [Nx, Ny, Nz] = GaussNormal(X,Y,Z,sigma)
    %GaussNormal computes smoothened normal using Gaussian derivative
    %
    % GaussNormal(X,Y,Z,sigma)
    
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

    [Xz, Xphi] = GaussGradient(X, sigma);
    [Yz, Yphi] = GaussGradient(Y, sigma);
    [Zz, Zphi] = GaussGradient(Z, sigma);

    Nx = Yz.*Zphi - Zz.*Yphi;
    Ny = Zz.*Xphi - Xz.*Zphi;
    Nz = Xz.*Yphi - Yz.*Xphi;
    Nnorm = sqrt(Nx.^2 + Ny.^2 + Nz.^2);
    Nx = Nx./Nnorm;
    Ny = Ny./Nnorm;
    Nz = Nz./Nnorm;

end

