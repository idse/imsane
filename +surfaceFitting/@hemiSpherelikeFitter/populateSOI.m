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
    %X0  = @(Z) polyval(this.fittedParam.pX, double(Z), this.fittedParam.SX, this.fittedParam.muX);
    %Y0  = @(Z) polyval(this.fittedParam.pY, double(Z), this.fittedParam.SY, this.fittedParam.muY);

    %-----------------------------------------------------------------
    % Deal with desired chart dependencies
    %-----------------------------------------------------------------
    
    desCharts = struct();
    for i = 1:length(this.charts)
        desCharts.(this.charts(i).name) = this.charts(i).desired;
    end
    
    
    
    %-----------------------------------------------------------------
    % Define cylinder1 
    %-----------------------------------------------------------------
        
    % define the image of the chart
    name       = 'cylinder';
    if desCharts.(name) == 1
        
        debugMsg(2, 'generating cylinder chart\n');
        
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
    % Calculate induced metric
    %-----------------------------------------------------------------

    % after having added the embeddings in different charts we can now
    % calculate the induced metric
    
    SOI.NCalcInducedMetric(gti);
	
    %-----------------------------------------------------------------
    % Proper charts cylinder1
    %-----------------------------------------------------------------
    
    origChartNames = {};
    
    if desCharts.('cylinder_proper') == 1 
        origChartNames = [origChartNames, 'cylinder'];
    end
 
    if ~isempty(origChartNames)
        
        for i = 1:length(origChartNames)
            
            debugMsg(2, ['generating ' origChartNames{i} ' chart\n']);

            origChart = SOI.atlas(gti).getChart(origChartNames{i});
            domain = origChart.domain;

            dz = origChart.image.stepSize(1);
            df = origChart.image.stepSize(2);
            
            domain.name

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
            
            
        end
    end
    
    
end
