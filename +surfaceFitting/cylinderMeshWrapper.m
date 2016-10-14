classdef cylinderMeshWrapper < surfaceFitting.meshWrapper
    % Generate a SurfaceOfInterest object with atlas containing cylinder 
    % charts from an externally produced triangular mesh representation of 
    % a surface. This is not possible for all surfaces.
       
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
    
    %---------------------------------------------------------------------
    % properties
    %--------------------------------------------------------------------- 
   
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = cylinderMeshWrapper()
            % MESHWRAPPER Create a SOI from a generic mesh
            %
            % meshWrapper()
            %

            
            % call superclass constructor
            
            this = this@surfaceFitting.meshWrapper();
            
            % overload charts, since second level child classes run into
            % issues with instantiating the abstract property charts; We
            % may want to not declare an abstract property in the
            % superclass surfaceFitter and therfore instead overload in all
            % children - This is a do to. 
            
            % charts - charts this fitter can produce and their properties
            %
            % the rows have the structure: 
            % {name, description, stepsize, intersections, fundamental, desired}
            %
            % stepsize:         2d vector of stepsize
            % intersections:    lists indices of charts with overlapping domain
            % fundamental:      boolean, does chart define a set in the topology
            % desired:          boolean indicating whether to produce it
            
            this.charts = struct(...
            'name', {'cylinder1', 'cylinder2',...
                        'cylinder1_proper', 'cylinder2_proper', ...
                        'equidistant'},...
            'description', {'Cylinder coordinates z, phi',...
                            'Cylinder shifted by pi',...
                            'Proper cylinder coordinates',...
                            'Shifted proper cylinder coordinates',...
                            'Equidistant maps from seedpoints, aka exponential map'},...
            'stepSize', {[1 .1], [1 .1], [1 1], [1 1], [1 1]},...
            'intersections', {[2 3 4],[1 3 4],[1 2],[1 2],[]},...%{[2 3 4], [1 3 4], [], [], [], []},...
            'fundamental', {1, 1, 0, 0, 1},...
            'desired', {1, 1, 0, 0, 0});
            
            % initialize fitOptions
            this.fitOptions.chartSeeds      = [];     
            this.fitOptions.transitionWidth = 100;
            this.fitOptions.fixAxis         = 1;
            this.fitOptions.rotation        = eye(4);
            this.fitOptions.translation     = eye(4);
            this.fitOptions.fixResolution   = 0;
            this.fitOptions.resolution      = [];
            
            % initialize fittedParam
            this.fittedParam = struct('mesh', [], 'submeshes', [],...
                                      'submVidx', [], 'maxDist',[],...
                                      'intersectIdx',{cell(0)},...
                                      'intersects',[]);
        end
        
        
        % ------------------------------------------------------
        %  generate embedding
        % ------------------------------------------------------
               
        function [embeddings, charts] = generateEmbedding(this, chartName)
            % GENERATEEMBEDDING Create the embedding in some chart
            %
            % [embedding, chart] = generateEmbedding(chartName)
            %
            % The embedding is generated using the result of a fit, stored
            % in fittedParam, so fitSurface should be called first.
            % The available charts are listed in the charts property of the
            % SurfaceFitter object.
            % 
            % See also surfaceFitting.tpsFitter.fitSurface
            
            % number of seeds
            %seeds = this.fitOptions.chartSeeds;
            %nSeeds = numel(seeds);
            
            desCharts = struct();
            for i = 1:length(this.charts)
                desCharts.(this.charts(i).name) = this.charts(i).desired;
            end
            
            embeddings = {};
            charts = {};
            
            % --------------------------------------------
            % Define Cylinder 1 & 2 Embedding
            % --------------------------------------------
            
            this.fittedParam.submeshes{1}.u = cell(1,sum(cat(1,this.charts.desired)));
            
            debugMsg(2, ['generateEmbedding(' chartName ')\n']);
            
            if strcmp(chartName, 'cylinder1') || strcmp(chartName, 'cylinder2') 
            
                % here, we parametrize the mesh using a cylindrical
                % coordinate system. For this we orient the data such that
                % the long axis coincides with z; 
                
                % extract the point cloud to determine its orientation. 
                
                points = this.fittedParam.mesh.v;
                
                if all([sum(sum((this.fitOptions.rotation-eye(4)).^2)),...
                        sum(sum((this.fitOptions.translation-eye(4)).^2))]==0) || ...
                        this.fitOptions.fixAxis == 0
                    
                    debugMsg(2, 'Updating orientation\n');
                    
                    pc = surfaceDetection.PointCloud(points);
                    pc.determineROI(5);
                    rotation    = pc.ROI.rotation;
                    translation = pc.ROI.translation;

                    this.fitOptions.rotation    = rotation;
                    this.fitOptions.translation = translation;
                end
                
                % from setRoi in the PointCloud class ;
                R = this.fitOptions.rotation;
                T = this.fitOptions.translation;
                
                pointsHomogenous = [points ones([size(points,1) 1])];
                alignedPoints = (R*T*pointsHomogenous')';
                alignedPoints = alignedPoints(:,1:3);
                
    
                z   = alignedPoints(:,3); %  u
                phi = atan2(alignedPoints(:,2),alignedPoints(:,1))+pi; % v
                r   = sqrt( (alignedPoints(:,1)).^2+(alignedPoints(:,2)).^2 );
                
                if strcmp(chartName, 'cylinder2')
                    phi = mod(phi + (phi<pi)*pi + (phi>pi)*pi,2*pi);
                end
                
                this.fittedParam.submeshes{1}.u{1} = [z(:),phi(:)];
                
                % determine step size in phi;
                r_max = max(r);
                %dz   = 1; % this will be the user fed input. 
                %dphi = dz/r_max; 
                
                c1step = this.getChartResolution(chartName);
                
                if isempty(this.fitOptions.resolution)
                    dz   = c1step(1);
                    dphi = dz/r_max;
                    this.fitOptions.resolution = [dz dphi];
                else
                    dz   = this.fitOptions.resolution(1);
                    dphi = this.fitOptions.resolution(2);
                end
                % define the stepsize in the charts; 
                stepSize = [dz dphi];
                this.setChartResolution(chartName, stepSize);
                
                % define the boundary of the future chart;
%                boundary = {[min(round(z)), max(round(z))], [0, 2*pi]};
                boundary = {[min(round(z)), max(round(z))], [min(phi) max(phi)]};
%                boundary = {[min(round(z)), max(round(z))],[pi/2 3*pi/2]};
                        
                [zG,phiG] = meshgrid(boundary{1}(1):dz:boundary{1}(2),...
                    boundary{2}(1):dphi:boundary{2}(2));
                % phiG is the phi coordinate in cylinder1; 
                x = alignedPoints(:,1);
                y = alignedPoints(:,2);
                   
                % finally we are in the position to generate the embedding
                % grids for x and y; 
                % we multiply everywhere with r_max to get phi into real
                % distances, rather than angles; Otherweise the matlab
                % internal routine to compute the triangulation produces
                % errors.
                if exist('scatteredInterpolant','file') > 0
                    F = scatteredInterpolant(z,r_max*phi,x,'linear','none');
                else 
                    F = TriScatteredInterp(z,r_max*phi,x);
                end
                xG = F(zG,r_max*phiG);

                if exist('scatteredInterpolant','file') > 0
                    F = scatteredInterpolant(z,r_max*phi,y,'linear','none');
                else 
                    F = TriScatteredInterp(z,r_max*phi,y);
                end
                yG = F(zG,r_max*phiG);
                
                % don't forget to rotate [xG,yG,zG] back into the embedding space!
            
                grid = [xG(:),yG(:),zG(:)];
                
                gridHomogenous = [grid ones([size(grid,1) 1])];
                alignedGrid = (inv(T)*inv(R)*gridHomogenous')';
                alignedGrid = alignedGrid(:,1:3);
                
                % this is the embeddding grid
                grids{1} = reshape(alignedGrid(:,1),size(xG,1),size(xG,2));
                grids{2} = reshape(alignedGrid(:,2),size(xG,1),size(xG,2));
                grids{3} = reshape(alignedGrid(:,3),size(xG,1),size(xG,2));
                
                
                % and now we deal with the ImSAnE infrastructure, i.e. make 
                % charts, embedding etc. ...
                
                % define the image and domain of the chart
                % boundary was defined above; so was stepsize
                chartim  = diffgeometry.FiniteSet(chartName,...
                                                boundary, stepSize);
                domain   = chartim.makeIndexSet();

                % chart: name_index -> name
                % the domain is the image of the chart, the image is the domain of
                % the fit, i.e. the embedding space
                chart = diffgeometry.CoordinateMap(domain, chartim,...
                                                chartim.makeHandles());

                % for charts made with FiniteSet.makeHandles we can also set an
                % analytic inverse (for speed and accuracy)
                chartInv = diffgeometry.CoordinateMap(chartim, domain,...
                                        chartim.makeInverseHandles());
                chart.setInverse(chartInv);

                % embedding: cylinder1 -> targetSpace
                embedding = diffgeometry.CoordinateMap(chartim,...
                                         this.fitDomain, grids);

                % embedding \circ chart: cylinder1_index -> targetSpace
                embedding = embedding.compose(chart); 

                %
                embeddings{1} = embedding;
                charts{1} = chart;
                
            end        
        end

        %------------------------------------------------------
        % populate SOI
        %------------------------------------------------------
        
        function  populateSOI(this, SOI, varargin)
            % POPULATESOI Add fit result to SurfaceOfInterest object
            %
            % populateSOI(SOI)
            % populateSOI(SOI, t)
            %
            % SOI:  SurfaceOfInterest object
            % t:    time, needs to be provided if SOI.dynamic = true
            %
            % Generates chart domain, chart and embedding using the result 
            % of a fit, adds these to the SOI.
            % Fit results are stored in fittedParam, so fitSurface should 
            % be called first.
            % 
            % See also surfaceFitting.tpsFitter.fitSurface
            
            % gti : geometric time index (one for static geometry, time index for
            % dynamic geometry)
            if SOI.dynamic == false
                gti = 1;
                if length(varargin) == 1
                    debugMsg(1, 'static SOI: ignoring time provided\n');
                end
            else
                if length(varargin) == 1
                    gti = SOI.tIdx(varargin{1});
                else
                    error('for dynamic SOI, time argument needs to be provided');
                end
            end
            
            desCharts = struct();
            for i = 1:length(this.charts)
                desCharts.(this.charts(i).name) = this.charts(i).desired;
            end

            %-----------------------------------------------------------------
            % Define cylinder1 
            %-----------------------------------------------------------------

            % WE ALWAYS GENERATE CYLINDER 1, OTHER CHARTS DEPEND ON IT
            
            name = 'cylinder1';
            debugMsg(2, ['generating ' name ' charts \n']);

            [embeddings, charts] = this.generateEmbedding(name);

%             % TODO: topology = overlap, transition maps
%             for i = 1:numel(embeddings)

            intersects = {};
            SOI.topologicalSpace(gti).addSet(charts{1}.domain, intersects);
            SOI.atlas(gti).addChart(charts{1});
            SOI.embedding(gti).addPatch(embeddings{1});
%            end

            %-----------------------------------------------------------------
            % Define cylinder2 
            %-----------------------------------------------------------------

            % define the image of the chart
            name       = 'cylinder2';
            if desCharts.(name) == 1

                debugMsg(2, 'generating cylinder2 chart\n');

                [embeddings, charts] = this.generateEmbedding(name);
                
                intersects = {'cylinder1_index'};
                SOI.topologicalSpace(gti).addSet(charts{1}.domain, intersects);
                SOI.atlas(gti).addChart(charts{1});
                SOI.embedding(gti).addPatch(embeddings{1});

                cyl1 = SOI.atlas(gti).getChart('cylinder1');
                cyl2 = SOI.atlas(gti).getChart('cylinder2');
                
                % generate transition map between first and second cylinder.    
                % analytic definition of the transition map: shift by - pi, mod 2 pi
                % mod is not actually needed the way it is implemented here
                phi_12   = { @(u)(u{1}), @(u)( u{2} + (u{2} < pi).*pi - (u{2} > pi).*pi ) };
                tmap12   = diffgeometry.CoordinateMap(cyl1.image, cyl2.image, phi_12);

                % generate transition map between second an first cylinder.    
                phi_21 = { @(u)(u{1}), @(u) mod( u{2} + (u{2} < pi).*pi + (u{2} > pi).*pi ,2*pi) };           
                tmap21   = diffgeometry.CoordinateMap(cyl2.image, cyl1.image, phi_21); 

                % add transition maps to the atlas
                SOI.atlas(gti).addTransitionMap(tmap12);
                SOI.atlas(gti).addTransitionMap(tmap21);
            end
            
            %-----------------------------------------------------------------
            % Define equidistant / exponential maps
            %-----------------------------------------------------------------

            if desCharts.('equidistant') == 1      

                nVorSeeds = numel(this.fitOptions.VorSeeds);
                nDisks = numel(this.fitOptions.diskSeeds);

                for j = nVorSeeds+1:nVorSeeds+nDisks

                    debugMsg(2, ['generating equidistant chart ' num2str(j) '\n']);

                    v = this.fittedParam.submeshes{j}.v;
                    f = this.fittedParam.submeshes{j}.f;

                    % align the mesh vertices
                    points = v;
                    R = this.fitOptions.rotation;
                    T = this.fitOptions.translation;
                    pointsHomogenous = [points ones([size(points,1) 1])];
                    alignedPoints = (R*T*pointsHomogenous')';
                    alignedPoints = alignedPoints(:,1:3);
                    v = alignedPoints;

    %                 % VISUALIZE:
    %                 v = xp.fitter.fittedParam.submeshes{1}.v;
    %                 f = xp.fitter.fittedParam.submeshes{1}.f;
    %                 fv.Vertices = [v(:,1), v(:,2), v(:,3)];
    %                 fv.Faces = f;
    %                 fv.FaceVertexCData = kron([1 1 1],mat2gray(v(:,3)));
    %                 patch(fv, 'EdgeColor','k', 'FaceColor', 'none');
    %                 shading flat;
    %                 axis equal;
    %                 view([0 0 -1]);

                    % map the full mesh seedIdx into the submesh
                    seedIdx = this.fitOptions.diskSeeds(j);
                    vidx = this.fittedParam.submVidx{j};
                    idxMap = double(0*vidx);
                    idxMap(vidx) = 1:sum(vidx);
                    seedIdx = idxMap(seedIdx);

                    % compute the normal at the seed
                    [normal,~] = compute_normal(v,f);
                    normal = normal';
                    seedNormal = normal(seedIdx,:);

                    % distance map
                    D = perform_fast_marching_mesh(v, f, seedIdx);

                    % seed normal is z direction
                    % construct x and y for phi
                    x = [1 0 0];
                    x = x - dot(seedNormal, x)*seedNormal;
                    x = x./norm(x);
                    y = cross(seedNormal, x);

                    % THIS MAY NEED SOME WORK
                    % I translate relative to seed but really I want distance
                    % to the axis, for a cylindrical sample its about the same
                    % though
                    nVerts = size(v,1);
                    vp = v - repmat(v(seedIdx,:),[nVerts 1]);
                    vx = dot(vp',repmat(x, [nVerts 1])')';
                    vy = dot(vp',repmat(y, [nVerts 1])')';
                    phi = atan2(vy, vx);

                    % back to cartesian
                    X = D.*cos(phi);
                    Y = D.*sin(phi);
                    u = [X Y];

                    % interpolate mesh coordinates over regular grid of surface 
                    % coordinates to obtain embedding
                    [uG, vG] = meshgrid(min(X):max(X), min(Y):max(Y));
                    embX = griddata(u(:,1), u(:,2), v(:,1), uG, vG);
                    embY = griddata(u(:,1), u(:,2), v(:,2), uG, vG);
                    embZ = griddata(u(:,1), u(:,2), v(:,3), uG, vG);

                    gridHomogenous = [embX(:) embY(:) embZ(:) ones([numel(embX) 1])];
                    alignedGrid = (inv(T)*inv(R)*gridHomogenous')';
                    alignedGrid = alignedGrid(:,1:3);

                    % this is the embeddding grid
                    embX = reshape(alignedGrid(:,1),size(embX,1),size(embX,2));
                    embY = reshape(alignedGrid(:,2),size(embX,1),size(embX,2));
                    embZ = reshape(alignedGrid(:,3),size(embX,1),size(embX,2));

                    % add everything to the SOI object

                    % index set
                    bdry = {[1, size(uG,2)], [1, size(uG,1)]};
                    domain = diffgeometry.FiniteSet(['equidistant_' num2str(j) '_index'], bdry, [1 1]);

                    % chart
                    bdry      = ({[min(uG(:)),max(uG(:))],[min(vG(:)),max(vG(:))]});
                    image     = diffgeometry.FiniteSet(['equidistant_' num2str(j)],bdry,[1 1]);
                    chart     = diffgeometry.CoordinateMap(domain, image, {uG, vG});

                    grids = {embX, embY, embZ};

                    embedding = diffgeometry.CoordinateMap(chart.domain, this.fitDomain,...
                                 grids); 

                    intersects = {};
                    SOI.topologicalSpace(gti).addSet(chart.domain, intersects);
                    SOI.atlas(gti).addChart(chart);
                    SOI.embedding(gti).addPatch(embedding);         
                end
            end

            %-----------------------------------------------------------------
            % Proper charts 
            %-----------------------------------------------------------------

            origChartNames = {};

            if desCharts.('cylinder1_proper') == 1 
                origChartNames = [origChartNames, 'cylinder1'];
            end

            if desCharts.('cylinder2_proper') == 1    
                origChartNames = [origChartNames, 'cylinder2'];
            end

            if ~isempty(origChartNames)
        
                % after having added the embeddings in different charts we can now
                % calculate the induced metric
                % here we need to 
                SOI.NCalcInducedMetric(SOI.timePoints(gti));

                for i = 1:length(origChartNames)

                    debugMsg(2, ['generating ' origChartNames{i} ' chart\n']);

                    origChart = SOI.atlas(gti).getChart(origChartNames{i});
                    domain = origChart.domain;

                    dz = origChart.image.stepSize(1);
                    df = origChart.image.stepSize(2);
                    
                    gzz = SOI.g(gti).getPatch(domain.name).cmp({1,1});
                    gff = SOI.g(gti).getPatch(domain.name).cmp({2,2});

                    % remove nans; 
                    gzz(isnan(gzz)) = 0;
                    gff(isnan(gff)) = 0;
                    
                    % zp = int_0^z dz sqrt(gzz) 
                    % z is first coordinate so second index
                    zp = real(cumsum(sqrt(gzz)*dz, 2));

                    [yy,xx] = meshgrid(1:size(zp,2),1:size(zp,1));
                    ind = ~isnan(zp);
                    if exist('scatteredInterpolant','file')>0
                        zpI = scatteredInterpolant(xx(ind),yy(ind),zp(ind),'linear','linear');
                    else
                        zpI = TriScatteredInterp(xx(ind),yy(ind),zp(ind));
                    end
                    zp = zpI(xx,yy);
                    
                    
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

                    this.fittedParam.submeshes{1}.u{2+i} = [zp(:),fp(:)];
                    
                    % use the same stepsize for the proper chart image
                    bdry = {[min(zp(:)), max(zp(:))], [min(fp(:)), max(fp(:))]};

                    image = diffgeometry.FiniteSet([origChartNames{i} '_proper'], bdry, [dz dz]);
                    pumpkin = diffgeometry.CoordinateMap(domain, image, {zp, fp});

                    SOI.atlas(gti).addChart(pumpkin);
                    
                end
            end
        end
        
    end
end