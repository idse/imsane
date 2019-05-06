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

    ux    = @(u) real( sqrt(rsq(u{1})).*cos(u{2}) );
    uy    = @(u) real( sqrt(rsq(u{1})).*sin(u{2}) );

    EmbeddingAnalytic = {@ (u) ux(u) + X0(u{1}), ...
            @(u) uy(u)+ Y0(u{1}), @(u) real(u{1})};
   
    % --------------------------------------------
    % Define Cylinder Embedding
	% --------------------------------------------
        
    name       = 'cylinder';
    if strcmp(chartName, name)
        
        debugMsg(2, ['generating ' name ' embedding\n']);
        
        % define the image of the chart
        boundary   = {[this.fittedParam.zmin, this.fittedParam.zmax], [this.fittedParam.phimin, this.fittedParam.phimax]};
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
end
