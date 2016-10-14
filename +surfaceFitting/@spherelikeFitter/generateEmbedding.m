function [embedding, chart] = generateEmbedding(this, chartName)
    % generate embedding for spherelikeFitter using parameterization of
    % surface. 
    %
    % [embedding, chart] = generateEmbedding(chartName);
    %
    % We fit sperelike surfaces in a cylindrical basis as described in the
    % supporting information. 
    
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
    
    % available charts are: Cylinder1, Cylinder2, PolarLowerZ, PolarUpperZ
    
   
    % functional form of the embedding is used for all embeddings
    
    % Embedding is an analytic map. Handles are:
    X0  = @(Z) polyval(this.fittedParam.pX, double(Z), this.fittedParam.SX, this.fittedParam.muX);
    Y0  = @(Z) polyval(this.fittedParam.pY, double(Z), this.fittedParam.SY, this.fittedParam.muY);
    rsq = @(Z) polyval(this.fittedParam.pR, double(Z), this.fittedParam.SR, this.fittedParam.muR);
    e   = @(Z) polyval(this.fittedParam.pe, double(Z), this.fittedParam.Se, this.fittedParam.mue); 

    
    alpha = this.fittedParam.medianphase;
    phase = @(Z) polyval(this.fittedParam.pphase,double(Z), this.fittedParam.Sphase,this.fittedParam.muphase)+alpha; 
    % A polygonal fit to the raw phase.
    %sinphase = @(Z) fitFourierEval(this,double(Z));
    % A fourier fit to the sine of the phase.
    %sinphase = @(Z) polyval(this.fittedParam.psinphase,double(Z), this.fittedParam.Ssinphase,this.fittedParam.musinphase); 
    % A polygonal fit to the sine of the phase.


    ux    = @(u) real( sqrt(rsq(u{1})).*cos(u{2}) );
    uy    = @(u) real(1./sqrt(1-e(u{1}).^2).*sqrt(rsq(u{1})).*sin(u{2}));

    EmbeddingAnalytic = {@ (u) cos(phase(u{1})).*ux(u) - sin(phase(u{1})).*uy(u) + X0(u{1}), ...
            @(u) sin(phase(u{1})).*ux(u) + cos(phase(u{1})).*uy(u)+ Y0(u{1}), @(u) real(u{1})};
    % Evalulation of the embedding when raw phase is fitted.    
    %EmbeddingAnalytic = {@ (u) real(sqrt(1-sinphase(u{1}).^2)).*ux(u) - sinphase(u{1}).*uy(u) + X0(u{1}), ...
    %       @(u) sinphase(u{1}).*ux(u) + real(sqrt(1-sinphase(u{1}).^2)).*uy(u)+ Y0(u{1}), @(u) real(u{1})};
    % Evalulation of the embedding when sine of phase is fitted. 
    
    % --------------------------------------------
    % Define Cylinder1 Embedding
	% --------------------------------------------
        
    name       = 'cylinder1';
    if strcmp(chartName, name)
        
        debugMsg(2, ['generating ' name ' embedding\n']);
        
        % define the image of the chart
        boundary   = {[this.fittedParam.zmin, this.fittedParam.zmax], [0, 2*pi]};
        stepSize   = this.getChartResolution(name);
        image      = diffgeometry.FiniteSet(name, boundary, stepSize);
        
        domain    = image.makeIndexSet();
        
        % chart: cylinder1_index -> cylinder1
        chart = diffgeometry.CoordinateMap(domain, image, image.makeHandles());
        % for charts made with FiniteSet.makeHandles we can also set an
        % analytic inverse (for speed and accuracy)
        chartInv = diffgeometry.CoordinateMap(image, domain, image.makeInverseHandles());
        chart.setInverse(chartInv);
        
        % embedding: cylinder1 -> targetSpace
        % the domain is the image of the chart, the image is the domain of
        % the fit, i.e. the embedding space
        embedding = diffgeometry.CoordinateMap(image,...
                                                    this.fitDomain, EmbeddingAnalytic);
                                                
        % embedding \circ chart: cylinder1_index -> targetSpace
        embedding = embedding.compose(chart); 
        
        return;
    end
    
    % --------------------------------------------
    % Define Cylinder2 embedding
    % --------------------------------------------
    
    name       = 'cylinder2';
    if strcmp(chartName, name)
                
        debugMsg(2, 'generating cylinder2 embedding\n');

        % cylinder2 embedding is defined as cylinder1 embedding shifted by
        % pi and is most easily obtained by a transition map composed with
        % cylinder1 embedding
        [~, chart1] = this.generateEmbedding('cylinder1');
                
        % since the second chart is just one with the cut in a different place the
        % coordinates can actually be the same: the transition map will specify
        % where the zeros are relative to each other
        
        boundary   = {[this.fittedParam.zmin, this.fittedParam.zmax], [0, 2*pi]};
        stepSize   = this.getChartResolution(name);
        image      = diffgeometry.FiniteSet(name, boundary, stepSize);     

        domain    = image.makeIndexSet();
        
        % chart: cylinder2_index -> cylinder2
        chart = diffgeometry.CoordinateMap(domain, image, image.makeHandles());
        chartInv = diffgeometry.CoordinateMap(image, domain, image.makeInverseHandles());
        chart.setInverse(chartInv);
        
        % generate transition map between second and first cylinder.    
        % analytic definition of the transition map: shift by pi, mod 2 pi
        phi_21 = { @(u)(u{1}), @(u) mod( u{2} + (u{2} < pi).*pi + (u{2} > pi).*pi ,2*pi) };           
        
        % tmap21 : cylinder2 -> cylinder1
        tmap21   = diffgeometry.CoordinateMap(chart.image, chart1.image, phi_21); 

        % compose transition map with chart2 to get cylinder1 coordinates on
        % cylinder2_index domain
        % t21_c1 : cylinder2_index -> cylinder1
        t21_c1 = tmap21.compose(chart);

        % embedding: cylinder1 -> targetSpace
        embedding = diffgeometry.CoordinateMap(chart1.image, this.fitDomain, EmbeddingAnalytic);
        
        % embedding : cylinder2_index -> targetSpace
        embedding = embedding.compose(t21_c1);
        
        return;
    end
    
    % --------------------------------------------
    % Define polarLowerZ embedding
    % -------------------------------------------- 
    
    name       = 'polarLowerZ';
    if strcmp(chartName, name)
        
        debugMsg(2, 'generating polarLowerZ embedding\n');
        
        % get the grids of the embedding of cylinder1
        embedding1Grids = this.generateEmbedding('cylinder1').apply();        
        
        % get z single valued over the xy plane
        ind = embedding1Grids{3} < this.fittedParam.zRmaxLower;
        X   = double(embedding1Grids{1}(ind));
        Y   = double(embedding1Grids{2}(ind));

        % Now determine xmin, ymin, xmax, ymax
        xmin = floor(min(X));
        xmax = ceil(max(X));
        ymin = floor(min(Y));
        ymax = ceil(max(Y));

        % define the image of the chart
        boundary   = {[xmin, xmax], [ymin, ymax]};
        stepSize   = this.getChartResolution(name);
        image      = diffgeometry.FiniteSet(name, boundary, stepSize); 
        
        domain    = image.makeIndexSet();
        
        % chart: polarLowerZ_index ->  polarLowerZ
        chart = diffgeometry.CoordinateMap(domain, image, image.makeHandles());
        chartInv = diffgeometry.CoordinateMap(image, domain, image.makeInverseHandles());
        chart.setInverse(chartInv);
        
        % we are using x,y as coordinates in the polar chart, which means
        % that z(x,y) is the only non-trivial component
        
        PLZGrids = chart.apply; % x,y
        PLZGrids{1} = double(PLZGrids{1});
        PLZGrids{2} = double(PLZGrids{2});

        zPLZ = getZOfXY(this.fittedParam.zmin, this.fittedParam.zRmaxLower);

        embedding = diffgeometry.CoordinateMap(chart.domain, this.fitDomain,...
            {PLZGrids{1},PLZGrids{2}, zPLZ});
        
        return;
    end
    
    % --------------------------------------------
    % Define polarUpperZ
    % --------------------------------------------
    
    name       = 'polarUpperZ';
    if strcmp(chartName, name)
        
        debugMsg(2, 'generating polarUpperZ embedding\n');
        
        % get the grids of the embedding of cylinder1
        embedding1Grids = this.generateEmbedding('cylinder1').apply();        
        
        % get z single valued over the xy plane
        ind = embedding1Grids{3} > this.fittedParam.zRmaxUpper;
        X   = double(embedding1Grids{1}(ind));
        Y   = double(embedding1Grids{2}(ind));

        % Now determine xmin, ymin, xmax, ymax
        xmin = floor(min(X));
        xmax = ceil(max(X));
        ymin = floor(min(Y));
        ymax = ceil(max(Y));

        % define the image of the chart
        boundary   = {[xmin, xmax], [ymin, ymax]};
        stepSize   = this.getChartResolution(name);
        image      = diffgeometry.FiniteSet(name, boundary, stepSize); 
        
        domain    = image.makeIndexSet();

        % chart: polarLowerZ_index ->  polarUpperZ
        chart = diffgeometry.CoordinateMap(domain, image, image.makeHandles());
        chartInv = diffgeometry.CoordinateMap(image, domain, image.makeInverseHandles());
        chart.setInverse(chartInv);
        
        % we are using x,y as coordinates in the polar chart, which means
        % that z(x,y) is the only non-trivial component
        
        PUZGrids = chart.apply; % x,y
        PUZGrids{1} = double(PUZGrids{1});
        PUZGrids{2} = double(PUZGrids{2});

        zPUZ = getZOfXY(this.fittedParam.zRmaxUpper, this.fittedParam.zmax);

        embedding = diffgeometry.CoordinateMap(chart.domain, this.fitDomain,...
            {PUZGrids{1}, PUZGrids{2}, zPUZ});
        
        return;
    end
            
	%--------------------------------------------
    % get z(x,y)
    %--------------------------------------------
    
    function zOfXY = getZOfXY(zmin, zmax)
        % a fast way to obtain Z(X,Y)
        %
        % zOfXY = getZOfXY(zmin, zmax)
        %
        % we have x = R(z) cos(phi) + x0(z) etc, and we want z(x,y)
        % we can first obtain z(R), then R(x,y)
        % R(x,y) is the solution of R^2 = (x-x0(R))^2 + (y-y0(R))^2*(1-e^2), 
        % we can rewrite that as x(R, y) = sqrt( R^2 - (y-y0(R))^2*(1-e^2) ) + x0(R) 
        % and for each y invert to get R(x,y) 

        % step 1: get z(R)
        zIndex  = zmin:zmax;
        rOfZ    = sqrt(rsq(zIndex));
        rIndex  = 0:floor(max(rOfZ));

        [rOfZ, idx, ~] = unique(rOfZ);        
        zIndex = zIndex(idx);
        zOfR    = interp1(rOfZ, zIndex, rIndex);	
        
        eofZOfR = e(zOfR);
        
        % step 2: x(R,y)
        xIndex = floor(this.fitDomain.boundary{1}(1)):ceil(this.fitDomain.boundary{1}(2));
        yIndex = floor(this.fitDomain.boundary{2}(1)):ceil(this.fitDomain.boundary{2}(2));
        xSize  = length(xIndex);
        ySize  = length(yIndex);
        
        X0OfR = X0(zOfR);
        Y0OfR = Y0(zOfR);
        
        rOfXY = zeros([ySize xSize]);
        zOfXY = NaN*ones([ySize xSize]);
        
        % ymin and ymax could potentially be negative. Using yindex for the
        % loop instead, as yy is used as an index for a matrix (zOfXY and 
        % others) later;

        for yy = 1:length(yIndex)%ymin:ymax
            
            yi = yIndex(yy);
            
            inroot = (rIndex.^2 - (yi - Y0OfR).^2.*(1-eofZOfR.^2));
            
            rootPve = inroot > 0;
            rPve = rIndex(rootPve);

            if length(rPve) > 2
                
                root = sqrt(inroot(rootPve));

                % for fixed y past together the two roots to get r(x,y)
                root1 = -root + X0OfR(rootPve);
                root2 = root + X0OfR(rootPve);

                % use spline interpolation to deal with the gap in the
                % middle
                % remove duplicates first so interpolation works
                    
                xOfRY = [fliplr(root1) root2];
                rPveBoth =  [fliplr(rPve) rPve];
                
                mask = diff(xOfRY) > 0;
                
                [xOfRY, idx, ~] = unique(xOfRY);
                rPveBoth =  rPveBoth(idx);

                rOfXY(yy,:) = interp1(xOfRY(mask), rPveBoth(mask), xIndex, 'spline');

                % mask to get rid of the spline crazyness outside the data
                xMask = xIndex < floor(max(xOfRY)) & xIndex > ceil(min(xOfRY));
                rOfXY(yy,~xMask) = NaN;

                % zOfXY also needs spline interpolation because it goes
                % beyond the zmax for which there is data for z(R)
                zOfXY(yy,:) = interp1(rOfZ, zIndex, rOfXY(yy,:), 'spline');
                zOfXY(yy,~xMask)  = NaN;
               
            end
        end
        
        % a few pixels of spline crazyness could have survived on the boundary
        while(max(abs(zOfXY(:))) > 2*max(abs(zmax),abs(zmin))) 
            mask = imerode(abs(zOfXY) > min(abs(zmax),abs(zmin)), strel('disk', 1));
            zOfXY(~mask) = NaN;
        end
        
	% Also remove values far below zmin;
	if sign(zmin) == 1
		while(min(abs(zOfXY(:))) < zmin/4) 
		    mask = imerode(abs(zOfXY) > min(abs(zmax),abs(zmin)), strel('disk', 1));
		    zOfXY(~mask) = NaN;
		    %disp('removing spline crazyness');
        end
%     else 
% 	     while(min(zOfXY(:)) < zmin) 
%             	mask = imerode(abs(zOfXY) > min(abs(zmax),abs(zmin)), strel('disk', 1));
%             	zOfXY(~mask) = NaN;
% 	        %disp('removing spline crazyness');
%          end	
	end
        %zOfXY = zOfXY(ymin:ymax, xmin:xmax);        
        ybd = [find(yIndex==ymin),find(yIndex==ymax)];
        xbd = [find(xIndex==xmin),find(xIndex==xmax)];
        zOfXY = zOfXY(ybd(1):ybd(2), xbd(1):xbd(2));
        
        
    end
end
