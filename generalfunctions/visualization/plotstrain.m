function plotstrain(s, xbar, ybar, options)
    % plot a strain structure
    %
    % plotstrain(s, xbar, ybar, options)
    %
    % s is the strain structure 
    % with fields MajorAxisLength, MinorAxisLength, Rotation
    % xbar, ybar are the coordinates of the position to draw the straincross
    % options: 'Rotation' 0,1
    % strainscale: some rescaling, default times 50

    
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
    
    
    theta = s.Orientation*pi/180;
    if isfield(options,'strainscale')
        scale = options.strainscale;
    else
        scale = 50;
    end

    N = 10;
    majAx = scale*abs(s.MajorAxisLength)/2;
    majort = -majAx:majAx/N:majAx;
    majorx = cos(theta)*majort;
    majory = -sin(theta)*majort;

    minAx = scale*abs(s.MinorAxisLength)/2;
    minort = -minAx:minAx/N:minAx;
    minorx = sin(theta)*minort;
    minory = cos(theta)*minort;

    lw = 2;
    poscol = [1 0 0];
    negcol = [0 1 1];

    if s.MajorAxisLength > 0
        plot(xbar + majorx, ybar + majory, 'Color', poscol,'LineWidth',lw);
    else
        plot(xbar + majorx, ybar + majory, 'Color', negcol,'LineWidth',lw);
    end
    if s.MinorAxisLength > 0
        plot(xbar + minorx, ybar + minory, 'Color', poscol,'LineWidth',lw);
    else
        plot(xbar + minorx, ybar + minory, 'Color', negcol,'LineWidth',lw);
    end

    % plot rotation
    if isfield(options,'Rotation') && options.Rotation
        rotscale=5;
        phi = [linspace(0,-rotscale*s.Rotation,N) linspace(pi,pi-rotscale*s.Rotation,N)];
        cosphi = cos(phi);
        sinphi = sin(phi);
        R = [ cos(theta)   sin(theta)
             -sin(theta)   cos(theta)];
        xy = [majAx*cosphi; majAx*sinphi];
        xy = R*xy;
        x = xy(1,:) + xbar;
        y = xy(2,:) + ybar;
        plot(x(1:N),y(1:N),'g','LineWidth',2);
        plot(x(N+1:end),y(N+1:end),'g','LineWidth',2);
    end

end