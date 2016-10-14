function  populateSOI(this, SOI, varargin)
    % populate SOI for spherelikeFitter
    %
    % populateSOI(SOI)
    % populateSOI(SOI, t)
    %
    % if SOI is dynamic, time needs to be provided
    
    % populate SOI object by computing desired charts and arrange them in
    % the strucutres from the diffgeometry section such as the atlas. 
    
    % gti : geometric time index (one for static geometry, time index for
    % dynamic geometry)
    
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
    
    if SOI.dynamic == false
        gti = 1;
        if length(varargin) == 1
            disp('static SOI: ignoring time provided');
        end
    else
        if length(varargin) == 1
            gti = SOI.tIdx(varargin{1});
        else
            error('for dynamic SOI, time argument needs to be provided');
        end
    end
    
    % part of the embedding parametrization that is required here:
    % polyeval can have trouble with single precision variables -> Convert
    % to double;
    X0  = @(Z) polyval(this.fittedParam.pX, double(Z), this.fittedParam.SX, this.fittedParam.muX);
    Y0  = @(Z) polyval(this.fittedParam.pY, double(Z), this.fittedParam.SY, this.fittedParam.muY);
    e   = @(Z) polyval(this.fittedParam.pe, double(Z), this.fittedParam.Se, this.fittedParam.mue); 

    %-----------------------------------------------------------------
    % Deal with desired chart dependencies
    %-----------------------------------------------------------------
    
    desCharts = struct();
    for i = 1:length(this.charts)
        desCharts.(this.charts(i).name) = this.charts(i).desired;
    end
    
    % cylinder2 requires cylinder1
    if desCharts.cylinder2 == 1
        desCharts.cylinder1 = 1;
    end
    
    % if polar is desired, both cylinders will be required
    if desCharts.polarLowerZ == 1 || desCharts.polarUpperZ == 1
        desCharts.cylinder1 = 1;
        desCharts.cylinder2 = 1;
    end
    
    %-----------------------------------------------------------------
    % Define cylinder1 
    %-----------------------------------------------------------------
        
    % define the image of the chart
    name       = 'cylinder1';
    if desCharts.(name) == 1
        
        debugMsg(2, 'generating cylinder1 chart\n');
        
        [embedding1, chart1] = this.generateEmbedding(name);
        
        intersects = {};
        SOI.topologicalSpace(gti).addSet(chart1.domain, intersects);
        SOI.atlas(gti).addChart(chart1);
        SOI.embedding(gti).addPatch(embedding1);        
        
        % we also set fitter.embedding to cylinder1 (for quality inspection
        % etc)
        this.fittedPoints = embedding1.apply;
    end
    
    %-----------------------------------------------------------------
    % Define cylinder2 
    %-----------------------------------------------------------------

    % define the image of the chart
    name       = 'cylinder2';
    if desCharts.(name) == 1
        
        debugMsg(2, 'generating cylinder2 chart\n');
        
        [embedding2, chart2] = this.generateEmbedding(name);
        
        intersects = {'cylinder1_index'};
        SOI.topologicalSpace(gti).addSet(chart2.domain, intersects);
        SOI.atlas(gti).addChart(chart2);
        SOI.embedding(gti).addPatch(embedding2);
        
        % generate transition map between first and second cylinder.    
        % analytic definition of the transition map: shift by - pi, mod 2 pi
        % mod is not actually needed the way it is implemented here
        phi_12   = { @(u)(u{1}), @(u)( u{2} + (u{2} < pi).*pi - (u{2} > pi).*pi ) };
        tmap12   = diffgeometry.CoordinateMap(chart1.image, chart2.image, phi_12);

        % generate transition map between second an first cylinder.    
        phi_21 = { @(u)(u{1}), @(u) mod( u{2} + (u{2} < pi).*pi + (u{2} > pi).*pi ,2*pi) };           
        tmap21   = diffgeometry.CoordinateMap(chart2.image, chart1.image, phi_21); 
        
        % add transition maps to the atlas
        SOI.atlas(gti).addTransitionMap(tmap12);
        SOI.atlas(gti).addTransitionMap(tmap21);
    end
    
    %-----------------------------------------------------------------
    % Define polarLowerZ 
    %-----------------------------------------------------------------
    
    name       = 'polarLowerZ';
    if desCharts.(name) == 1
        
        debugMsg(2, 'generating polarLowerZ chart\n');
        
        [embeddingPLZ, chartPLZ] = this.generateEmbedding(name);
       
        intersects = {'cylinder1_index', 'cylinder2_index'};
        SOI.topologicalSpace(gti).addSet(chartPLZ.domain, intersects);
        SOI.atlas(gti).addChart(chartPLZ);
        SOI.embedding(gti).addPatch(embeddingPLZ);

        % generate transition map between polarLowerZ and cylinder1    
        % this takes x,y to z, phi
        % z(x,y) is already part of the embedding, phi(x,y) still has to be
        % computed
        
        PLZGrids = embeddingPLZ.apply;

        % the embedding is defined as
        % X = R cos(phi) + X0, Y = R sin(phi)/sqrt(1-e^2) + Y0
        cosphi = PLZGrids{1} - X0(PLZGrids{3});
        sinphi = sqrt( 1- e(PLZGrids{3}).^2 ).*( PLZGrids{2} - Y0(PLZGrids{3}) );
        % this phi is the angle, the next phi is a transition map!
        phiPLZ = atan2(sinphi, cosphi);
        
        % tmapLZ1 : polarLowerZ (x,y) -> cylinder1 (z, phi)
        phi_PLZ1 = { PLZGrids{3}, phiPLZ };           
        tmapPLZ1   = diffgeometry.CoordinateMap(chartPLZ.image, chart1.image, phi_PLZ1); 

        % now generate the transition map from cylinder1 to polarLowerZ,
        % taking z,phi -> x,y, this is trivial because it is just the
        % cylinder1 embedding over the lower half
        
        embedding1Grids = embedding1.apply;
        
        x   = embedding1Grids{1};
        x(embedding1Grids{3}> this.fittedParam.zRmaxLower) = NaN;
        y   = embedding1Grids{2};
        y(embedding1Grids{3}> this.fittedParam.zRmaxLower) = NaN;
        
        % tmapLZ1 : cylinder1 -> polarLowerZ
        phi_1PLZ = {x,y};
        tmap1PLZ   = diffgeometry.CoordinateMap(chart1.image, chartPLZ.image, phi_1PLZ); 
        
        % add both transition maps to the atlas
        SOI.atlas(gti).addTransitionMap(tmapPLZ1);
        SOI.atlas(gti).addTransitionMap(tmap1PLZ);

        % mow define transition map between cylinder2 and polarLowerZ 
        % this is much shorter, since we can use the tmap from cylinder1
        % to cylinder2. 

        % tmapLZ2: polarLowerZ -> cylinder2
        tmapPLZ2   = tmap12.compose(tmapPLZ1);
        % tmap2LZ: cylinder2 -> polarLowerZ
        tmap2PLZ   = tmap1PLZ.compose(tmap21);
        
        % add to the atlas
        SOI.atlas(gti).addTransitionMap(tmapPLZ2);
        SOI.atlas(gti).addTransitionMap(tmap2PLZ);
    end
    
    %-----------------------------------------------------------------
    % Define polarUpperZ
    %-----------------------------------------------------------------

	name       = 'polarUpperZ';
    if desCharts.(name) == 1
        
        debugMsg(2, 'generating polarUpperZ chart\n');
        
        [embeddingPUZ, chartPUZ] = this.generateEmbedding(name);
       
        intersects = {'cylinder1_index', 'cylinder2_index'};
        SOI.topologicalSpace(gti).addSet(chartPUZ.domain, intersects);
        SOI.atlas(gti).addChart(chartPUZ);
        SOI.embedding(gti).addPatch(embeddingPUZ);

        % generate transition map between polarUpperZ and cylinder1    
        % this takes x,y to z, phi
        % z(x,y) is already part of the embedding, phi(x,y) still has to be
        % computed
        
        PUZGrids = embeddingPUZ.apply;

        % the embedding is defined as
        % X = R cos(phi) + X0, Y = R sin(phi)/sqrt(1-e^2) + Y0
        cosphi = real(PUZGrids{1} - X0(PUZGrids{3}));
        sinphi = real(sqrt( 1- e(PUZGrids{3}).^2 ).*( PUZGrids{2} - Y0(PUZGrids{3}) ));
        % this phi is the angle, the next phi is a transition map!
        phiPUZ = atan2(sinphi, cosphi);
        
        % tmapLZ1 : polarUpperZ (x,y) -> cylinder1 (z, phi)
        phi_PUZ1 = { PUZGrids{3}, phiPUZ };           
        tmapPUZ1   = diffgeometry.CoordinateMap(chartPUZ.image, chart1.image, phi_PUZ1); 

        % now generate the transition map from cylinder1 to polarUpperZ,
        % taking z,phi -> x,y, this is trivial because it is just the
        % cylinder1 embedding over the upper half
        
        embedding1Grids = embedding1.apply;
        
        x   = embedding1Grids{1};
        x(embedding1Grids{3} < this.fittedParam.zRmaxUpper) = NaN;
        y   = embedding1Grids{2};
        y(embedding1Grids{3} < this.fittedParam.zRmaxUpper) = NaN;
        
        % tmapLZ1 : cylinder1 -> polarLowerZ
        phi_1PUZ = {x,y};
        tmap1PUZ   = diffgeometry.CoordinateMap(chart1.image,...
                                chartPUZ.image, phi_1PUZ); 
        
        % add both transition maps to the atlas
        SOI.atlas(gti).addTransitionMap(tmapPUZ1);
        SOI.atlas(gti).addTransitionMap(tmap1PUZ);

        % mow define transition map between cylinder2 and polarLowerZ 
        % this is much shorter, since we can use the tmap from cylinder1
        % to cylinder2. 

        % tmapLZ2: polarLowerZ -> cylinder2
        tmapPUZ2   = tmap12.compose(tmapPUZ1);
        % tmap2LZ: cylinder2 -> polarLowerZ
        tmap2PUZ   = tmap1PUZ.compose(tmap21);
        
        % add to the atlas
        SOI.atlas(gti).addTransitionMap(tmapPUZ2);
        SOI.atlas(gti).addTransitionMap(tmap2PUZ);
    end
    
	%-----------------------------------------------------------------
    % Calculate induced metric
    %-----------------------------------------------------------------

    % after having added the embeddings in different charts we can now
    % calculate the induced metric
    
    SOI.NCalcInducedMetric(gti);
	
    %-----------------------------------------------------------------
    % Proper charts cylinder1
    %-----------------------------------------------------------------
    
    origChartNames = {};
    
    if desCharts.('cylinder1_proper') == 1 
        origChartNames = [origChartNames, 'cylinder1'];
    end
    
    if desCharts.('cylinder2_proper') == 1    
        origChartNames = [origChartNames, 'cylinder2'];
    end
    
    if ~isempty(origChartNames)
        
        for i = 1:length(origChartNames)
            
            debugMsg(2, ['generating ' origChartNames{i} ' chart\n']);

            origChart = SOI.atlas(gti).getChart(origChartNames{i});
            domain = origChart.domain;

            dz = origChart.image.stepSize(1);
            df = origChart.image.stepSize(2);

            gzz = SOI.g(gti).getPatch(domain.name).cmp({1,1});
            gff = SOI.g(gti).getPatch(domain.name).cmp({2,2});

            % zp = int_0^z dz sqrt(gzz) 
            % z is first coordinate so second index
            zp = real(cumsum(sqrt(gzz)*dz, 2)); 
            
            % non-monotonous behaviour of z' as a function of z for phip ==
            % const lead to a re-definiton of z': z' = int
            % sqrt(gzz(z,phi_with_max_distance)) dz;
            % Identify the phip with maximum distance zp;
            [~,maxind] =  max(zp(:));
            [II,~] = ind2sub(size(zp),maxind);
            % generate a grid that replicates this zp; 
            [zp,~] = meshgrid(zp(II,:),1:size(zp,1));

            % fp = int_0^f df sqrt(gff) + c(z) (see egggeometry)
            intdphi = real(cumsum(sqrt(gff)*df, 1));
            mid = round(size(gff,1)/2); % mid ~ pi
            cz = pi - real(sum(sqrt(gff(1:mid,:))*df, 1));
            cz = repmat(cz, size(intdphi,1), 1);
            fp = intdphi + cz; 

            % use the same stepsize for the proper chart image
            bdry = {[min(zp(:)), max(zp(:))], [min(fp(:)), max(fp(:))]};

            image = diffgeometry.FiniteSet([origChartNames{i} '_proper'], bdry, [dz dz]);
            pumpkin = diffgeometry.CoordinateMap(domain, image, {zp, fp});

            SOI.atlas(gti).addChart(pumpkin);
            
            %-----------------------------------------------------------------
            % Proper charts at poles
            %-----------------------------------------------------------------
            % These charts are based on the proper distances we get through
            % the integration of the metric during generation of 
            % cylinder1_proper charts. To avoid recomputation, we put the
            % generation of anteriorEquidistant and posteriorEquidistant
            % charts in this position. 
            if i == 1 
                
                if desCharts.('anteriorEquidistant') == 1
            
                %%%%%
                % anteriorEquidistant charts;         
                %%%%%                        
                    
                debugMsg(2, 'generating anteriorEquidistant chart\n');
                
              
                % define a new chart, that contains the z' as z and phi as phi; 
                phi      = origChart.apply{2};
                z        = origChart.apply{1};
                eGrids   = SOI.embedding.getPatch(origChart.domain.name).apply();

                
                
                % perform a tri scattered interpolation to obtain z(zp,phi)
                % and so on.
        

                
                z_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(z(:)));
                x_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{1}(:)));
                y_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{2}(:)));
                
                
                % now define a meshgrid xG,yG, compute rG and phiG (raduis and angle)
                % in this meshgrid and plug thisinto z_zp_phi, thus relating
                % zp with r_G
                
                fac = .6;% defines reach of charts in fraction of max(zp). .6 is a good value 
                % for overlap with other charts.
                
                [xG,yG] = meshgrid(-max(zp(:))*fac:dz:max(zp(:))*fac);
                % get r and phi
                [phiG,rG] = cart2pol(xG,yG);
                phiG = phiG+pi;
                phiG(phiG>max(phi(:))) = max(phi(:));
                
                zG =  z_zp_phi(double(rG),double(phiG(end:-1:1,end:-1:1)));
                xGE = x_zp_phi(double(rG),double(phiG(end:-1:1,end:-1:1)));
                yGE = y_zp_phi(double(rG),double(phiG(end:-1:1,end:-1:1)));
                %zG  = z_zp_phi(double(rG),double(phiG));
                %xGE = x_zp_phi(double(rG),double(phiG));
                %yGE = y_zp_phi(double(rG),double(phiG));
                 
                clear z_zp_phi;
                clear x_zp_phi;
                clear y_yp_phi;
                % this is the definition of the embedding.
                zG  = griddata(double(xG(~isnan(zG))),double(yG(~isnan(zG))),double(zG(~isnan(zG))),double(xG),double(yG));
                xGE = griddata(double(xG(~isnan(xGE))),double(yG(~isnan(xGE))),double(xGE(~isnan(xGE))),double(xG),double(yG));
                yGE = griddata(double(xG(~isnan(yGE))),double(yG(~isnan(yGE))),double(yGE(~isnan(yGE))),double(xG),double(yG));
                
                zG(isnan(zG)) = 0;
                
                bdry = {[1, size(zG,1)], [1, size(zG,2)]};
                domain = diffgeometry.FiniteSet('anteriorEquidistant_index', bdry, [1 1]);
            
                bdry = {[-max(zp(:))*fac, max(zp(:))*fac], [-max(zp(:))*fac, max(zp(:))*fac]};
                image = diffgeometry.FiniteSet('anteriorEquidistant', bdry, [1 1]);
                chartPLZ = diffgeometry.CoordinateMap(domain, image, {xG, yG});
                
                bdry      = {[1, size(zG,1)], [1, size(zG,2)]};
                domain    = diffgeometry.FiniteSet('anteriorEquidistant_index', bdry, [1 1]);

                bdry      = ({[min(xG(:)),max(xG(:))],[min(yG(:)),max(yG(:))]});
                image     = diffgeometry.FiniteSet('anteriorEquidistant',bdry,[1 1]);
                chart     = diffgeometry.CoordinateMap(domain, image, {xG, yG});


                embeddingPLZ = diffgeometry.CoordinateMap(chart.domain, this.fitDomain,...
                             {xGE, yGE, zG});
                         
                intersects = {};
                SOI.topologicalSpace(gti).addSet(chartPLZ.domain, intersects);
                SOI.atlas(gti).addChart(chartPLZ);
                SOI.embedding(gti).addPatch(embeddingPLZ);         
                         
                
                % generate transition map between polarLowerZ2 and cylinder1    
                % this takes x,y to z, phi
                % z(x,y) and, phi(x,y) still have to be
                % computed

                % the embedding is defined as
                % X = R cos(phi) + X0, Y = R sin(phi)/sqrt(1-e^2) + Y0
                cosphi = xGE - X0(zG);
                sinphi = sqrt( 1- e(zG).^2 ).*( yGE - Y0(zG) );
                % this phi is the angle, the next phi is a transition map!
                phiPLZ = atan2(sinphi, cosphi);
                
                % there seems to be an offset between this phi and the phi
                % of the cylinder1 chart. We solve this by determination of
                % a phi in cylinder1, and finding the best match in the
                % two embeddings. The difference to phiPLZ will be
                % considered an offset. 
                bla = chart1.apply();
                ff = bla{2}(round(end/2),100);
                
                embedding1Grids = embedding1.apply();
                
                diff = (xGE-embedding1Grids{1}(round(end/2),100)).^2+(yGE-embedding1Grids{2}(round(end/2),100)).^2+...
                        (zG-embedding1Grids{3}(round(end/2),100)).^2;
                [~,ii] = min(diff(:));
                
                phiPLZ = mod(phiPLZ,2*pi);
                phioff = ff - phiPLZ(ii);
                
                phiPLZ = mod(phiPLZ+phioff,2*pi);

                % tmapLZ1 : polarLowerZ (x,y) -> cylinder1 (z, phi)
                phi_PLZ1 = { zG, phiPLZ};           
                tmapPLZ1   = diffgeometry.CoordinateMap(chartPLZ.image, chart1.image, phi_PLZ1); 

                % now generate the transition map from cylinder1 to polarLowerZ2,
                % By composing pumkin with inverse cylinder, we get a map from z into zp, 
                % then the relations of the x,y plane are simply zp*cos(phi). 
              
                bla = origChart.apply();
                tmap11P = pumpkin.composeInverse(chart1);
                bla2 = tmap11P.apply();
                x = (bla2{1}).*cos(bla{2});
                y = (bla2{1}).*sin(bla{2});
                

                % tmapLZ1 : cylinder1 -> polarLowerZ
                phi_1PLZ = {x,y};
                tmap1PPLZ   = diffgeometry.CoordinateMap(origChart.image, chartPLZ.image, phi_1PLZ); 
                
                % now obtain the transition map between cylinder and
                % pumpkin. Then compose the result with the transition map
                % from pumpkin to PLZ to get a transition map from cylinder
                % to PLZ.
                
                % curChart is the pumpkin1 Chart, chart1 the cylinder1
                %tmap11P = pumpkin.composeInverse(chart1);
                %tmap1PLZ = tmap1PPLZ.compose(tmap11P);
                
                
                % add both transition maps to the atlas
                SOI.atlas(gti).addTransitionMap(tmapPLZ1);
                SOI.atlas(gti).addTransitionMap(tmap1PPLZ);

                % mow define transition map between cylinder2 and polarLowerZ2 
                % this is much shorter, since we can use the tmap from cylinder1
                % to cylinder2. 

                % tmapLZ2: polarLowerZ -> cylinder2
                tmapPLZ2   = tmap12.compose(tmapPLZ1);
                % tmap2LZ: cylinder2 -> polarLowerZ
                tmap2PLZ   = tmap1PPLZ.compose(tmap21);

                % add to the atlas
                SOI.atlas(gti).addTransitionMap(tmapPLZ2);
                SOI.atlas(gti).addTransitionMap(tmap2PLZ);
                
                end

                
                %%%%%
                % posteriorEquidistant charts;         
                %%%%% 
                
                
                if desCharts.('posteriorEquidistant') == 1      
                         
                debugMsg(2, 'generating posteriorEquidistant chart\n');
                    
                phi      = origChart.apply{2};
                z        = origChart.apply{1};
                eGrids   = SOI.embedding.getPatch(origChart.domain.name).apply();
                
                zp = zp(:,end:-1:1);
                % we now formally define a chart with an inverted proper z.
                % this won't be added to the atlas. 
                bdry = {[min(zp(:)), max(zp(:))], [min(fp(:)), max(fp(:))]};

                image = diffgeometry.FiniteSet([origChartNames{i} '_proper_inv'], bdry, [dz dz]);
                pumpkin_semiInv = diffgeometry.CoordinateMap(origChart.domain, image, {zp, fp});
                
                
                z_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(z(:)));
                x_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{1}(:)));
                y_zp_phi = TriScatteredInterp(double(zp(:)),double(phi(:)),double(eGrids{2}(:)));
                
                
                fac = .6;% defines reach of charts in fraction of max(zp). .6 is a good value 

                % now define a meshgrid xG,yG, compute rG and phiG in this meshgrid and plug this 
                % into z_zp_phi
                [xG,yG] = meshgrid(-max(zp(:))*fac:dz:max(zp(:))*fac);
                % get r and phi
                [phiG,rG] = cart2pol(xG,yG);
                phiG = phiG+pi;
                phiG(phiG>max(phi(:))) = max(phi(:));
                
                zG =  z_zp_phi(double(rG),double(phiG(end:-1:1,end:-1:1)));
                xGE = x_zp_phi(double(rG),double(phiG(end:-1:1,end:-1:1)));
                yGE = y_zp_phi(double(rG),double(phiG(end:-1:1,end:-1:1)));
                %zG =  z_zp_phi(double(rG),double(phiG));
                %xGE = x_zp_phi(double(rG),double(phiG));
                %yGE = y_zp_phi(double(rG),double(phiG));
                 
                clear z_zp_phi;
                clear x_zp_phi;
                clear y_yp_phi;
                zG  = griddata(double(xG(~isnan(zG))),double(yG(~isnan(zG))),double(zG(~isnan(zG))),double(xG),double(yG));
                xGE = griddata(double(xG(~isnan(xGE))),double(yG(~isnan(xGE))),double(xGE(~isnan(xGE))),double(xG),double(yG));
                yGE = griddata(double(xG(~isnan(yGE))),double(yG(~isnan(yGE))),double(yGE(~isnan(yGE))),double(xG),double(yG));
                
                zG(isnan(zG)) = 0;
                
                bdry = {[1, size(zG,1)], [1, size(zG,2)]};
                domain = diffgeometry.FiniteSet('posteriorEquidistant_index', bdry, [1 1]);
            
                bdry = {[-max(zp(:))*fac, max(zp(:))*fac], [-max(zp(:))*fac, max(zp(:))*fac]};
                image = diffgeometry.FiniteSet('posteriorEquidistant', bdry, [dz dz]);
                chartPUZ = diffgeometry.CoordinateMap(domain, image, {xG, yG});
                
                bdry      = {[1, size(zG,1)], [1, size(zG,2)]};
                domain    = diffgeometry.FiniteSet('posteriorEquidistant_index', bdry, [1 1]);

                bdry      = ({[min(xG(:)),max(xG(:))],[min(yG(:)),max(yG(:))]});
                image     = diffgeometry.FiniteSet('posteriorEquidistant',bdry,[1 1]);
                chart     = diffgeometry.CoordinateMap(domain, image, {xG, yG});


                embeddingPUZ = diffgeometry.CoordinateMap(chart.domain, this.fitDomain,...
                             {xGE, yGE, zG});   
                
                         
                intersects = {};
                SOI.topologicalSpace(gti).addSet(chartPUZ.domain, intersects);
                SOI.atlas(gti).addChart(chartPUZ);
                SOI.embedding(gti).addPatch(embeddingPUZ);         
                         
      
                
                % Transition maps
                cosphi = xGE - X0(zG);
                sinphi = sqrt( 1- e(zG).^2 ).*( yGE - Y0(zG) );
                % this phi is the angle, the next phi is a transition map!
                phiPUZ = atan2(sinphi, cosphi);
                phiPUZ = mod(phiPUZ,2*pi);
                phi_PUZ1 = { zG, phiPUZ};           
                tmapPUZ1   = diffgeometry.CoordinateMap(chartPUZ.image, chart1.image, phi_PUZ1); 
                         
                 % now generate the transition map from cylinder1 to polarLowerZ2,
                % By composing pumkin with inverse cylinder, we get a map from z into zp, 
                % then the relations of the x,y plane are simply zp*cos(phi). 
              
                bla = origChart.apply();
                tmap11P = pumpkin_semiInv.composeInverse(chart1);
                bla2 = tmap11P.apply();
                x = (bla2{1}).*cos(bla{2});
                y = (bla2{1}).*sin(bla{2});
                
               
                
                
                % tmapUZ1 : cylinder1 -> polarUpperZ
                phi_1PUZ = {x,y};
                tmap1PPUZ   = diffgeometry.CoordinateMap(origChart.image, chartPUZ.image, phi_1PUZ); 
 
                SOI.atlas(gti).addTransitionMap(tmapPUZ1);
                SOI.atlas(gti).addTransitionMap(tmap1PPUZ);
                
                % mow define transition map between cylinder2 and polarUpperZ2 
                % this is much shorter, since we can use the tmap from cylinder1
                % to cylinder2. 

                % tmapUZ2: polarUpperZ -> cylinder2
                tmapPUZ2   = tmap12.compose(tmapPUZ1);
                % tmap2UZ: cylinder2 -> polarLowerZ
                tmap2PUZ   = tmap1PPUZ.compose(tmap21);

                % add to the atlas
                SOI.atlas(gti).addTransitionMap(tmapPUZ2);
                SOI.atlas(gti).addTransitionMap(tmap2PUZ);
                end
            end
        end
    end
    
    
end
