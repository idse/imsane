function plotMetricStrain(SOI, chartName, gridSize, strainScale)
    % plot the metric strain for a chart over the pullback in the first
    % channel to visualize map distortion
    %
    % plotMetricStrain(SOI, chartName, gridSize, strainScale)
    %
    % SOI:          SurfaceOfInterest object
    % chartName:    name of chart 
    % gridSize:     number of cells for which to calculate the strain
    % strainScale:  scale for visualization

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

    chart = SOI.atlas.getChart(chartName);
    chartDomain = chart.domain.name;
    Ngridlines = gridSize;
    
    gpatch = SOI.g.getPatch(chartDomain);
    if gpatch == 0
        SOI.NCalcInducedMetric;
        gpatch = SOI.g.getPatch(chartDomain);
    end
    
    g = gpatch.getTransform(chartName);

    chartGrids = chart.apply;
    
    u1 = chartGrids{1} - chart.image.boundary{1}(1);
    u2 = chartGrids{2} - chart.image.boundary{2}(1);

    du1 = chart.image.stepSize(1);
    du2 = chart.image.stepSize(2);

    di = (size(u1,1)-1)/Ngridlines(1);
    dj = (size(u1,2)-1)/Ngridlines(2);
    
    pb = SOI.data.getPatch(chartDomain).getTransform(chartName).apply{1};
    plotstrainOptions = struct('strainscale', strainScale);
    
    imshow(pb,[],'InitialMagnification',66)
    hold on 
    for i = 1:Ngridlines(1)
        for j = 1:Ngridlines(2)

            xi = round(i*di - di/2);
            yi = round(j*dj - dj/2);

            x = round(u1(xi,yi)/du1) + 1;
            y = round(u2(xi,yi)/du2) + 1;

            % the metric has to be transformed by the stepsize
            gGrids = g.apply;
            M = [   gGrids{1,1}(y, x)*du1^2    gGrids{1,2}(y, x)*du2*du1;
                    gGrids{2,1}(y, x)*du2*du1  gGrids{2,2}(y, x)*du2^2; ];

            if ~any(isnan(M(:))) && ~any(abs(M(:)) > 20) 

                s = metricStrainFromG(M);
                plotstrain(s, x, y, plotstrainOptions);
            end
        end
    end
    hold off
end