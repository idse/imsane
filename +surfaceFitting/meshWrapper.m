classdef meshWrapper < surfaceFitting.surfaceFitter
    % Generate a SurfaceOfInterest object with atlas containing conformal 
    % charts from an externally produced triangular mesh representation of 
    % a surface.
       
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
    
    properties (SetAccess = protected)
        
        % charts - charts this fitter can produce and their properties
        %
        % the rows have the structure: 
        % {name, description, stepsize, intersections, fundamental, desired}
        %
        % stepsize:         2d vector of stepsize
        % intersections:    lists indices of charts with overlapping domain
        % fundamental:      boolean, does chart define a set in the topology
        % desired:          boolean indicating whether to produce it
        charts = struct(...
            'name', {'exponential', 'conformal'},...
            'description', {'the exponential map', 'conformal map'},...
            'stepSize', {[1 1], [1 1]},...
            'intersections', {[]},...
            'fundamental', {1, 1},...
            'desired', {0, 1});
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = meshWrapper()
            % MESHWRAPPER Create a SOI from a generic mesh
            %
            % meshWrapper()
            %
            % vertices: Nx3 vertex position
            % faces:    Mx3 triangles
            
            % call superclass constructor
            this = this@surfaceFitting.surfaceFitter();
            
            % initialize fitOptions
            this.fitOptions.VorSeeds = [];     
            this.fitOptions.transitionWidth = 100;
            this.fitOptions.diskSeeds = [];
            this.fitOptions.diskRadius = 100;
            this.fitOptions.makeTMaps = true;
            this.fitOptions.fixedPtX = [];
            this.fitOptions.fixedPtU = [];
            
            % initialize fittedParam
            this.fittedParam = struct(  'mesh', [], 'submeshes', [],...
                                        'submVidx', [], 'submFidx', [],...
                                        'maxDist', [],...
                                        'intersectIdx',[], 'intersects',[]);
        end
        
        % ------------------------------------------------------
        % fitting
        % ------------------------------------------------------
        
        function fitSurface(this, mesh)
            % FITSURFACE Fit the pointcloud to produce an embedding patch
            % for tpsFitter using a thin plate spline fit
            %
            % fitSurface(mesh)
            %
            % mesh: struct with fields:
            %       v:    Nx3 vertex positions
            %       vn:   Nx3 vertex normals
            %       f:    Mx3 faces (triangles)
            %
            % fitOptions can be set through through the generic function
            % setFitOptions, the meshWrapper defines the following:
            %
            % VorSeeds:         seeds (vertex index) for the geodesic 
            %                   Voronoi cells
            % transitionWidth:  half-width of overlap between Voronoi cells
            % diskSeeds:        center (vertex index) for disk like submeshes
            % diskRadius:       radius of disk like submeshes
            % makeTMaps:        boolean, should transition maps be made?
            %
            % see also surfaceFitting.surfaceFitter.setFitOptions

            v = mesh.v;
            f = mesh.f;
            vn = mesh.vn;
            
            % make sure normals are normalized
            normalnorm = sqrt(vn(:,1).^2 + vn(:,2).^2 + vn(:,3).^2);
            for i=1:3
                vn(:,i) = vn(:,i)./normalnorm;
            end
            mesh.vn = vn;
            
            % store the mesh in fittedParam
            this.fittedParam.mesh = mesh;
            
            % fittedPoints are mesh vertices for now
            this.fittedPoints = {mesh.v(:,1), mesh.v(:,2), mesh.v(:,3)};
            
            %--------------------------------------------------------------
            % generate submeshes for patches on which charts will be defined
            %--------------------------------------------------------------

            seeds = this.fitOptions.VorSeeds;
            overlap = this.fitOptions.transitionWidth;
            nseeds = length(seeds);

            % distance map from each chart seed
            D = {};
            for i = 1:nseeds
                D{i} = perform_fast_marching_mesh(v, f, seeds(i,1))';
            end

            % vertex indices of submesh are those for which the distance to
            % a seedpoint are less than those to another plus the overlap
            submVidx = {};
            for i = 1:nseeds
                chartIdx = true(size(D{1}));
                for j = setdiff(1:nseeds,i)
                    chartIdx = chartIdx & D{i} <= D{j} + overlap;
                end
                submVidx{i} = chartIdx;
                
                % store the maximal distance between seed and edge of patch
                this.fittedParam.maxDist(i) = max(D{i}(submVidx{i}));
            end
            
            % disk submeshes
            diskSeeds = this.fitOptions.diskSeeds;
            nDiskSeeds = size(diskSeeds, 1);
            for i = 1:nDiskSeeds
                
                dist = perform_fast_marching_mesh(v, f, diskSeeds(i,1))';
                chartIdx = true(size(dist));
                for j = 1:nDiskSeeds
                    chartIdx = chartIdx & dist < this.fitOptions.diskRadius;
                end
                submVidx{nseeds + i} = chartIdx;
                
                % store the maximal distance between seed and edge of patch
                this.fittedParam.maxDist(nseeds + i) = max(dist(chartIdx));
            end

            % store in fittedParam
            this.fittedParam.submVidx = submVidx;

            % make the submeshes and store in fittedParam
            for i = 1:(nseeds+nDiskSeeds)
                
                [subm, vidx, fidx]= submeshFromV(mesh, submVidx{i});
                this.fittedParam.submVidx{i} = vidx;
                this.fittedParam.submFidx{i} = fidx;
                
                % store the seed with submesh indexing for later reference
                allSeeds = [seeds(:)' diskSeeds(:)'];
                subm.seed = sum(submVidx{i}(1:allSeeds(i)));
                
                this.fittedParam.submeshes{i} = subm;
            end
            
            % if no seeds make a trivial submesh 
            if nseeds + nDiskSeeds == 0
                this.fittedParam.submeshes{1} = mesh;
                this.fittedParam.submVidx{1} = 1:size(v,1); 
            end

            %--------------------------------------------------------------
            % determine overlap between submeshes
            %--------------------------------------------------------------
            
            intersects = false(nseeds+nDiskSeeds);
            intersectIdx = {};
            for i = 1:(nseeds+nDiskSeeds)
                
                vidxi = this.fittedParam.submVidx{i};
                
                for j = i+1:(nseeds+nDiskSeeds)
                    
                    vidxj = this.fittedParam.submVidx{j};
                    intersectIdx{i,j} = vidxi & vidxj;

                    % for practical reasons (interpolation to a smooth 
                    % transition map) we will want a respectable overlap
                    % when we call to sets intersecting, say, 100 triangles
                    if sum(intersectIdx{i,j}) > 100
                        intersects(i,j) = true;
                        intersects(j,i) = true;
                    end
                end
            end
            this.fittedParam.intersectIdx = intersectIdx;
            this.fittedParam.intersects = intersects;
        end
        
        %------------------------------------------------------
        % generate detector seed
        %------------------------------------------------------
        
        function seed = generateSeed(this)
            % GENERATESEED
            %
            % generateSeed()
            %
            % returns the mesh as a seed for the next surface detection
            
            seed = this.fittedParam.mesh;
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
            
            % seeds
            seeds = this.fitOptions.VorSeeds;
            nVorSeeds = size(this.fitOptions.VorSeeds,1);
            seeds = cat(1, seeds, this.fitOptions.diskSeeds);
            nCharts = numel(this.fittedParam.submeshes);
            
            embeddings = {};
            charts = {};
            
            for i = 1:nCharts
                                
                maxdist = this.fittedParam.maxDist(i);
                subm = this.fittedParam.submeshes{i};
                
                % seeds refer to complete mesh we may need the seed
                % index in the submeshes
                submVidx = this.fittedParam.submVidx{i};
                old2new = zeros(size(submVidx));
                old2new(submVidx) = 1:sum(submVidx);
                submSeed = old2new(seeds(i,1));
                
                %--------------------------------------------------
                % conformal maps
                %--------------------------------------------------
                
                if strcmp(chartName, 'conformal')

                    debugMsg(2, ['generating ' chartName ' embedding ' num2str(i) '\n']);
                    
                    % convert to structure used by matlabmesh toolbox
                    %subm_tb = makeMesh(subm.v, subm.f, subm.vn);
                    %u = embedSCP(subm_tb, 'robust'); %'fiedler', 'generalized'
                    
                    if isempty(this.fitOptions.fixedPtX)
                        
                        u = conformalParametrization(subm);
                        ucenter = u(submSeed,:);
                        u(:,1) = u(:,1) - ucenter(1);
                        u(:,2) = u(:,2) - ucenter(2);
                        
                    else
                        fixedPtX = this.fitOptions.fixedPtX{i};
                        fixedPtU = this.fitOptions.fixedPtU{i};
                        
                        fixedPtIdx = pointMatch(fixedPtX, subm.v);
                        constrVidx = [subm.seed fixedPtIdx];
                        constrU = [0 0; fixedPtU];

                        u = conformalParametrization(subm, constrVidx, constrU);
                    end

                    % rescale range of conformal map
                    scale = 2*maxdist/(max(u(:)) - min(u(:)));
                    u = u*scale;

                    % store vertex coordinates on submesh
                    subm.u{1} = u;

                %--------------------------------------------------
                % exponential maps
                %--------------------------------------------------
                
                elseif strcmp(chartName, 'exponential')
                    
                    debugMsg(2, ['generating ' chartName ' embedding ' num2str(i) '\n']);
                    [subm.u{1}, ~] = expmap(subm.v, subm.vn, submSeed);
                else
                    
                    error('nonexistent chart name');
                end
                
                % store vertex coordinates on submesh
                this.fittedParam.submeshes{i}.u{1} = u;

                %--------------------------------------------------
                % convert to grids, create imsane objects
                %--------------------------------------------------
                
                % Interpolate on grid
                u = subm.u{1}(:,1);
                v = subm.u{1}(:,2);

                notInf = ~isinf(u) & ~isinf(v);
                u = u(notInf);
                v = v(notInf);

                du = 1;
                dv = 1;
                
                minu = round(min(u(:)));
                maxu = round(max(u(:)));
                minv = round(min(v(:)));
                maxv = round(max(v(:)));

                [ugrid, vgrid] = meshgrid(minu:du:maxu, minv:dv:maxv);

                uidx = 1;
                mask = mesh_mask(subm, uidx, minv, minu, maxv-minv+1, maxu-minu+1);
                
                grids = {};
                for embIdx = 1:3

                    X = subm.v(notInf, embIdx);
                    if exist('scatteredInterpolant') > 0
                        F = scatteredInterpolant(u,v,X,'linear','none');
                    else
                        F = TriScatteredInterp(u,v,X);
                    end
                    grid = F(ugrid, vgrid);
                    grid(~mask) = NaN;
                    grids{embIdx} = grid;
                end

                % define the image and domain of the chart
                boundary = {[minu, maxu], [minv, maxv]};
                stepSize = [du dv];
                chartim  = diffgeometry.FiniteSet([chartName '_' num2str(i)],...
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

                embeddings{i} = embedding;
                charts{i} = chart;
            end
            
%             % store the rotation seeds for later reference
%             this.fitOptions.VorSeeds(:,2) = seeds(1:nVorSeeds, 2);
%             this.fitOptions.diskSeeds(:,2) = seeds(nVorSeeds+1:end, 2);
        end
        
        % ------------------------------------------------------
        %  populate SOI
        % ------------------------------------------------------
        
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

            %----------------------------------------------------------
            % Define charts
            %----------------------------------------------------------

            ncharts = numel(this.fittedParam.submeshes);
            
            % define the image of the chart
            for j = 1:numel(this.charts)
                
                name = this.charts(j).name;
                
                if desCharts.(name) == 1

                    debugMsg(2, ['generating ' name ' charts \n']);
                    [embeddings, charts] = this.generateEmbedding(name);
                    
                    for i = 1:ncharts
                        
                        % when adding to atlas, we can only specify
                        % intersections with sets that were added
                        % previously
                        chartIdx = find(this.fittedParam.intersects(:,i));
                        chartIdx = chartIdx(chartIdx < i);
                        isectNames = {};
                        for c = chartIdx'
                           isectNames = [isectNames, [name '_' num2str(c) '_index']];
                        end
                        
                        SOI.topologicalSpace(gti).addSet(charts{i}.domain, isectNames);
                        SOI.atlas(gti).addChart(charts{i});
                        SOI.embedding(gti).addPatch(embeddings{i});
                        
                        %-----------------------
                        % transition maps
                        %-----------------------
                        
                        if this.fitOptions.makeTMaps

                            disp('generating transition maps');

                            for c = chartIdx'

                                disp(['between ' num2str(c) ' and ' num2str(i)]);

                                % order of c, i matters, see makeTMaps
                                [tmapij, tmapji] = this.makeTMaps(c, i, charts{c}, charts{i});
                                SOI.atlas(gti).addTransitionMap(tmapij);
                                SOI.atlas(gti).addTransitionMap(tmapji);
                            end
                        end
                    end
                end
            end
        end
        
        % ------------------------------------------------------
        %  Evolve surface along surface normal.
        % ------------------------------------------------------

        function normallyEvolve(this, shift)
            % NORMALLYEVOLVE Normally shift surface out or in by some amount
            %
            % normallyEvolve(shift)
            %
            % shift : pixel distance along the normal
            
            disp('bla');
            
            v = this.fittedParam.mesh.v;
            vn = this.fittedParam.mesh.vn;
            
            v = v + shift*vn;
            this.fittedParam.mesh.v = v;
            this.fittedPoints = {v(:,1), v(:,2), v(:,3)};
            
            % also update submesh vertex positions
            for i = 1:(length(this.fitOptions.VorSeeds) + length(this.fitOptions.diskSeeds))
                disp('normally evolving submeshes');
                this.fittedParam.submeshes{i}.v = v(this.fittedParam.submVidx{i},:);
            end
        end

        %------------------------------------------------------
        % smoothing
        %------------------------------------------------------
        
        function smoothMesh(this, laplacian_type, T)
            % SMOOTHMESH smooth the mesh using heat diffusion
            %
            % smoothMesh(laplacianType, smoothingStrength)
            %
            % laplacianType:    'conformal', 'distance', 'combinatorial'
            % smoothingStrengh: positive number
            %
            % calls methods from Gabriel Peyre's toolbox_graph

            options.Tmax = T;
            
            vertex = this.fittedParam.mesh.v';
            face = this.fittedParam.mesh.f;

            % Laplacian
            if not(strcmp(laplacian_type,'conformal'))
                options.symmetrize = 1;
                options.normalize = 0;
                L = compute_mesh_laplacian(vertex,face,laplacian_type,options);
            else
                options.symmetrize = 0;
                options.normalize = 1;
                L = compute_mesh_laplacian(vertex,face,laplacian_type,options);
            end
            
            options.dt = 0.01;
            smoothv = perform_mesh_heat_diffusion(vertex,face,L,options);
            smoothv = smoothv';
            
            % update mesh
            this.fittedParam.mesh.v = smoothv;
            this.fittedPoints = {smoothv(:,1), smoothv(:,2), smoothv(:,3)};
            
            % also update submesh vertex positions
            for i = 1:length(this.fittedParam.submeshes)
                this.fittedParam.submeshes{i}.v = smoothv(this.fittedParam.submVidx{i},:);
            end
        end
        
        %------------------------------------------------------
        % visualization
        %------------------------------------------------------
        
        function inspectMesh(this, submeshIndex)
            % INSPECTTPS 3D visualization of thin plate spline fit
            %
            % inspectMesh()
            % inspectMesh(submeshIndex)
            %
            % submeshIndex: integer selecting submesh
            
            v = this.fittedParam.mesh.v;
            n = length(this.fittedParam.submeshes);
            colors = hsv(n);
            seeds = this.fitOptions.VorSeeds;
            
            if nargin==1 
                fv.Vertices = v;
                fv.Faces = this.fittedParam.mesh.f;
                fv.FaceVertexCData = kron([1 1 1],mat2gray(v(:,3)));
            else
            
                if numel(submeshIndex) > 1
                    fv.Vertices = [v(:,1), v(:,2), v(:,3)];
                    fv.Faces = this.fittedParam.mesh.f;
                    fv.FaceVertexCData = zeros(size(v));
                    for i = submeshIndex
                        fv.FaceVertexCData = fv.FaceVertexCData + kron(colors(i,:),...
                            this.fittedParam.submVidx{i}); 
                    end
                else
                    subm = this.fittedParam.submeshes{submeshIndex};
                    fv.Vertices = [subm.v(:,1), subm.v(:,2), subm.v(:,3)];
                    fv.Faces = subm.f;
                    fv.FaceVertexCData = kron([1 1 1],mat2gray(subm.v(:,3)));
                end
            end
            
            clf;
            patch(fv);
            shading interp;
            camlight;
            axis equal;
            view([1 0 0]);
            
            if nargin == 1
                hold on
                scatter3(v(seeds,1),v(seeds,2),v(seeds,3),'.r', 'fill');
                hold off
            end
        end
    end
    
    %---------------------------------------------------------------------
    % private methods
    %---------------------------------------------------------------------
    
    methods (Access = protected)
        
        %------------------------------------------------------
        % transition maps
        %------------------------------------------------------
        
        function [tmapij, tmapji] = makeTMaps(this, i, j, charti, chartj)
            % generate transition map
            %
            % [tmapij, tmapji] = makeTMaps(i, j, charti, chartj)
            %
            % i, j:             seed indices, assumes i < j
            % charti, chartj:   corresponding charts
            % tmapij, tmapji:   transition maps
          
            % some bullshit bookkeeping to get logical indices into each submesh for
            % the intersection with the other
            if this.fittedParam.intersects(i,j)

                submVidxi = this.fittedParam.submVidx{i};
                old2newi = zeros(size(submVidxi));
                old2newi(submVidxi) = 1:sum(submVidxi);

                submVidxj = this.fittedParam.submVidx{j};
                old2newj = zeros(size(submVidxj));
                old2newj(submVidxj) = 1:sum(submVidxj);

                isectIdx = this.fittedParam.intersectIdx{i,j};

                isectIdxi = false([sum(submVidxi) 1]);
                isectIdxi(old2newi(isectIdx)) = true;

                isectIdxj = false([sum(submVidxj) 1]);
                isectIdxj(old2newj(isectIdx)) = true;
            end

            % now generate subsubmeshes for the intersection
            submi = this.fittedParam.submeshes{i};
            submj = this.fittedParam.submeshes{j};

            submij = submeshFromV(submi, isectIdxi, false); 
            submji = submeshFromV(submj, isectIdxj, false);
            
            % visualize as pointcloud
            %figure, scatter3(submij.u(:,1), submij.u(:,2), submji.u(:,1));
            %figure, plot(submij.u(submij.b,1), submij.u(submij.b,2));
            %figure, scatter(submij.u(submij.b,1), submij.u(submij.b,2));
            %axis equal;

            % get the domains for charts on the two submeshes, needed to make an ImSAnE
            chartImi = charti.image;
            chartImj = chartj.image;

            % exclude interpolation outside the mesh
            uidx = 1;
            maski = mesh_mask(submij, uidx, chartImi.boundary{2}(1),...
                             chartImi.boundary{1}(1),...
                             chartImi.gridSize(2), chartImi.gridSize(1));

            maskj = mesh_mask(submji, uidx, chartImj.boundary{2}(1),...
                             chartImj.boundary{1}(1),...
                             chartImj.gridSize(2), chartImj.gridSize(1));

            % interpolate over chart image grids
            ugridsi = chartImi.makeGrids;
            ugridsj = chartImj.makeGrids;

            joveri = {};
            ioverj = {};
            for k = 1:2

                uij = submij.u{1};
                uji = submji.u{1};
                
                if exist('scatteredInterpolant') > 0
                    F = scatteredInterpolant(uij(:,1), uij(:,2), uji(:,k), 'linear', 'none');
                else
                    F = TriScatteredInterp(uij(:,1), uij(:,2), uji(:,k));
                end
                joveri{k} = F(double(ugridsi{1}), double(ugridsi{2}));
                joveri{k}(~maski) = NaN;

                if exist('scatteredInterpolant') > 0
                    F = scatteredInterpolant(uji(:,1), uji(:,2), uij(:,k), 'linear', 'none');
                else
                    F = TriScatteredInterp(uji(:,1), uji(:,2), uij(:,k));
                end
                ioverj{k} = F(double(ugridsj{1}), double(ugridsj{2}));
                ioverj{k}(~maskj) = NaN;
            end

            tmapij   = diffgeometry.CoordinateMap(chartImi, chartImj, joveri);
            tmapji   = diffgeometry.CoordinateMap(chartImj, chartImi, ioverj);
        end
    end
end