function geom = GaussGeometry(X,Y,Z,sigma)
    % GaussGeometry computes geometry using Gaussian derivative
    %
    % geom = GaussGeometry(X,Y,Z,sigma)
    %
    % X,Y,Z: embedding grids
    % sigma: width of Gaussian derivative
    %
    % geom is structure with fields
    %
    % g:    metric
    % ginv: inverse metric
    % dA:   area element ( = sqrt(detg))
    % N:    normal
    % K:    Gaussian curvature
    % H:    extrinsic curvature

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
    
    X = {X,Y,Z};

    [Xu, Xv] = GaussGradient(X{1}, sigma);
    [Yu, Yv] = GaussGradient(X{2}, sigma);
    [Zu, Zv] = GaussGradient(X{3}, sigma);

    N = {};
    N{1} = Yu.*Zv - Zu.*Yv;
    N{2} = Zu.*Xv - Xu.*Zv;
    N{3} = Xu.*Yv - Yu.*Xv;
    
    g = {};
    g{1,1} = Xu.*Xu + Yu.*Yu + Zu.*Zu;
    g{2,1} = Xu.*Xv + Yu.*Yv + Zu.*Zv;
    g{1,2} = g{2,1};
    g{2,2} = Xv.*Xv + Yv.*Yv + Zv.*Zv;

    detg = g{1,1}.*g{2,2} - g{1,2}.^2;

    ginv = {};
    ginv{1,1} = g{2,2}./detg;
    ginv{2,2} = g{1,1}./detg;
    ginv{1,2} = -g{1,2}./detg;
    ginv{2,1} = ginv{1,2};

    dA = sqrt(N{1}.^2 + N{2}.^2 + N{3}.^2); % = sqrt(detg)
    
    for i = 1:3, N{i} = N{i}./dA; end
    
    [Xuu, Xvv] = GaussGradient(Xu, sigma);
    [Yuu, Yuv] = GaussGradient(Yu, sigma);
    [Zuu, Zuv] = GaussGradient(Zu, sigma);

    [~, Xvv] = GaussGradient(Xv, sigma);
    [~, Yvv] = GaussGradient(Yv, sigma);
    [~, Zvv] = GaussGradient(Zv, sigma);

    L = {};
    L{1,1} = Xuu.*N{1} + Yuu.*N{2} + Zuu.*N{3};
    L{1,2} = Xvv.*N{1} + Yuv.*N{2} + Zuv.*N{3};
    L{2,1} = L{1,2};
    L{2,2} = Xvv.*N{1} + Yvv.*N{2} + Zvv.*N{3};

    H = (ginv{1,1}.*L{1,1} + ginv{2,2}.*L{2,2})/2;
    K = L{1,1}.*L{2,2} - L{1,2}.^2;

    geom = struct('g',{g},'ginv',{ginv},'dA',dA,'N',{N},'H',H,'K',K);
end

