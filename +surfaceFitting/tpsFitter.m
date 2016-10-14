classdef tpsFitter < surfaceFitting.surfaceFitter
    % Generate surface from a point cloud using thin plate spline fit
    %
    % The thin plate spline method fits a surface that behaves like a plate
    % with some bending stiffness. This surface remains parametrized by
    % coordinates xy on a plane, which we are also calling planar, so this
    % method should be used to fit point clouds that correspond to planar
    % surfaces.
       
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
    
    properties (Access = protected)
        
        redPointCloud;  % reduced point cloud forming input for the fit
    end
    
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
            'name', {'xy'},...
            'description', {'over the xy plane'},...
            'stepSize', {[1 1]},...
            'intersections', {[]},...
            'fundamental', {1},...
            'desired', {1});
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = tpsFitter()
            % TPSFITTER Create a new thin plate spline fitter
            %
            % tpsFitter()
            
            % call superclass constructor
            this = this@surfaceFitting.surfaceFitter();
            
            % initialize fitOptions
            this.fitOptions.smoothing = 1000;
            this.fitOptions.gridSize = [50 50];
            this.fitOptions.fitMask = [];
            
            % initialize fittedParam
            this.fittedParam = struct('tpsGridX',[],'tpsGridY',[],'tpsGridZ',[]);
        end
        
        % ------------------------------------------------------
        % fitting
        % ------------------------------------------------------
        
        function fitSurface(this, pointCloud)
            % FITSURFACE Fit the pointcloud to produce an embedding patch
            % for tpsFitter using a thin plate spline fit
            %
            % fitSurface(pointCloud)
            %
            % pointCloud:   a PointCloud object
            %
            % fitOptions can be set through through the generic function
            % setFitOptions, the tpsFitter defines the following:
            %
            % gridSize:     size of grid on which to generate fitted surface
            %               default [50 50], full size takes long
            % smoothing:    TPS smoothing parameter (default 1000) 
            % fitmask:      throw out the the resulting fit outside this mask (e.g.
            %               outside the convex hull of the data, where spline goes
            %               crazy)
            %
            % see also surfaceFitting.surfaceFitter.setFitOptions

            % set fit domain
            bdry = {round(pointCloud.ROI.xpRange),...
                round(pointCloud.ROI.ypRange), round(pointCloud.ROI.zpRange)};
            this.fitDomain = diffgeometry.FiniteSet('alignedEmbeddingSpace',...
                                                            bdry, [1 1 1]);
                                                        
            % tps fit works well with ~10^3 points, not too much more
            ssexp = floor(log10(size(pointCloud.points,1))) - 3;
            PC = pointCloud.points(1:10^ssexp:end,:);
            
            debugMsg(2, 'Fitting thin plate spline');
            ticID = tic;
            
            this.TPSfit(PC);
            
            dt = toc(ticID);
            debugMsg(2, ['dt = ' num2str(dt) ' sec\n']);
            
            this.fittedPoints = this.generateEmbedding('xy').apply;
        end
        
        %------------------------------------------------------
        % generate detector seed
        %------------------------------------------------------
        
        function seed = generateSeed(this)
            % GENERATESEED
            %
            % generateSeed()
            %
            % returns the Z component of the embedding as a seed for the
            % next surface detection
            
            seed = this.fittedPoints{3};
        end
        
        % ------------------------------------------------------
        %  generate embedding
        % ------------------------------------------------------
               
        function [embedding, chart] = generateEmbedding(this, chartName)
            % GENERATEEMBEDDING Create the embedding in some chart
            %
            % generateEmbedding(chartName)
            %
            % The embedding is generated using the result of a fit, stored
            % in fittedParam, so fitSurface should be called first.
            % The available charts are listed in the charts property of the
            % SurfaceFitter object.
            % 
            % See also surfaceFitting.tpsFitter.fitSurface
            
            name  = 'xy';
            if strcmp(chartName, name)
                
                debugMsg(2, ['generating ' name ' embedding\n']);
                
                % define the image and domain of the chart
                boundary   = this.fitDomain.boundary(1:2);
                stepSize   = this.getChartResolution(name);
                image      = diffgeometry.FiniteSet(name, boundary, stepSize);
                domain    = image.makeIndexSet();
                
                % chart: cylinder1_index -> cylinder1
                chart = diffgeometry.CoordinateMap(domain, image, image.makeHandles());
                % for charts made with FiniteSet.makeHandles we can also set an
                % analytic inverse (for speed and accuracy)
                chartInv = diffgeometry.CoordinateMap(image, domain, image.makeInverseHandles());
                chart.setInverse(chartInv);

                % interpolate tpsfit on image sized grid
                XYgrids = image.makeGrids();
                X = XYgrids{1};
                Y = XYgrids{2};

                Z = interp2(this.fittedParam.tpsGridX, this.fittedParam.tpsGridY,...
                                        this.fittedParam.tpsGridZ, X, Y);
                                                
                % cut off the fit when it goes outside the data
                zmin = this.fitDomain.boundary{3}(1);
                zmax = this.fitDomain.boundary{3}(2);

                Z(Z > zmax) = NaN;
                Z(Z < zmin) = NaN;
                
                if ~isempty(this.fitOptions.fitMask)
                    Z(~this.fitOptions.fitMask) = NaN;
                end

                % 3-4-14: line below creates holes normally evolving up
                % Z(Z < this.fitDomain.boundary{3}(1)) = NaN;
                grids = {X, Y, Z}; 
                
                % embedding: cylinder1 -> targetSpace
                embedding = diffgeometry.CoordinateMap(image,...
                                         this.fitDomain, grids);

                % embedding \circ chart: cylinder1_index -> targetSpace
                embedding = embedding.compose(chart); 
            end
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
            % Define xy
            %----------------------------------------------------------

            % define the image of the chart
            name       = 'xy';
            if desCharts.(name) == 1

                debugMsg(2, ['generating ' name ' chart\n']);
                [embedding1, chart1] = this.generateEmbedding(name);

                % this is necessary if we want fitter.embedding to be defined
                % when loading from disc
                this.fittedPoints = embedding1.apply;
                
                intersects = {};
                SOI.topologicalSpace(gti).addSet(chart1.domain, intersects);
                SOI.atlas(gti).addChart(chart1);
                SOI.embedding(gti).addPatch(embedding1);
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
            
            % create a pointcloud from the surface
            X = {this.fittedParam.tpsGridX, this.fittedParam.tpsGridY, this.fittedParam.tpsGridZ};
            
            % compute surface normal
            [Nx,Ny,Nz] = surfnorm(X{1}, X{2}, X{3});
            
            % point cloud in 3 x np matrix
            x  = [X{1}(:), X{2}(:), X{3}(:)];
            % normal in 3 x np matrix
            n  = [Nx(:),Ny(:),Nz(:)];
            % shift point cloud by normal.
            xs = x + shift*n;

            % remove nans;
            for k = 1 : 3
                n(isnan(x(:,k)),:)   = 0;
                x(isnan(x(:,k)),:)   = 0;
                xs(isnan(xs(:,k)),:) = 0;
            end

            newPointCloud = double(xs);
            
            % generate the new surface from this
            Z = griddata(newPointCloud(:,1), newPointCloud(:,2),...
                newPointCloud(:,3), double(this.fittedParam.tpsGridX),...
                double(this.fittedParam.tpsGridY));
           
            this.fittedParam.tpsGridZ = Z;
            this.fittedPoints = this.generateEmbedding('xy').apply;
        end
        
        % ------------------------------------------------------
        %  Evolve surface along z
        % ------------------------------------------------------

        function zEvolve(this, shift)
            % ZEVOLVE Shift surface along the z direction
            %
            % zEvolve(shift)
            %
            % shift:    pixels by which to shift
            % 
            % This is similar to normallyEvolve, but often more relevant to
            % data stacks whose PSF is anisotropic along z.
            
            disp('bla');
            
            this.fittedParam.tpsGridZ = this.fittedParam.tpsGridZ + shift;
            this.fittedPoints = this.generateEmbedding('xy').apply;
        end
        
        % ------------------------------------------------------
        % visualization
        % ------------------------------------------------------
        
        function inspectTPS(this, varargin)
            % INSPECTTPS 3D visualization of thin plate spline fit
            %
            % inspectTPS()
            % inspectTPS(pointCloud)
            % 
            % pointCloud:   PointCloud object
            
            grids = this.fittedPoints;
            surf(double(grids{1}), double(grids{2}), double(grids{3}))
            colormap(jet);
            shading interp;
            axis equal;
            
            if nargin == 2
                
                PC = varargin{1}.points;
                hold on;
                plot3(PC(1:100:end,1), PC(1:100:end,2), PC(1:100:end,3), '*k')
                hold off;
            end
        end
    end
      
    %---------------------------------------------------------------------
    % private methods
    %---------------------------------------------------------------------
    
    methods (Access = private)
        
        function TPSfit(this, PC)
            % TPSFIT Fit surface using smoothing thin plate split
            %
            % TPSfit(PC)
            % 
            % PC is an 3xN array of point cloud values [x y z]
            
            lambda = this.fitOptions.smoothing; % smoothing parameter lambda
            
            % "potential" matrix between the data points
            rij = @(x,xp) sqrt(...
                (x(1,:)'*ones(1,length(xp))-ones(length(x),1)*xp(1,:)).^2 ...
                + (x(2,:)'*ones(1,length(xp))-ones(length(x),1)*xp(2,:)).^2);
            
            % generate kernel
            % the x,y values of the control points are in a 2xN matrix X
            % z values are in a matrix Z
            X = PC(:,1:2)';
            Z = PC(:,3)';
            r = rij(X,X);
            M = r.^2.*log(r);
            M(r == 0)=0;
            N = cat(2, ones(length(X),1), X');
           
            % solve for surface parameters
            Q = M + lambda*eye(length(X));
            NpQ = N'/Q;
            b = (NpQ*N)\(NpQ*Z');
            a = Q\(Z'-N*b);
            
            % grid on which to generate the surface
            % creating the distance matrix rxX between data points and
            % gridpoints on which to generate is costly, which is why we
            % generate some smaller grid and later interpolate that
            % linearly
            xmin = this.fitDomain.boundary{1}(1);
            xmax = this.fitDomain.boundary{1}(2);
            ymin = this.fitDomain.boundary{2}(1);
            ymax = this.fitDomain.boundary{2}(2);
            
            Ngrid = this.fitOptions.gridSize;
            [tpsGridX, tpsGridY] = meshgrid(...
                                linspace(xmin, xmax, Ngrid(1)),...
                                linspace(ymin, ymax, Ngrid(2)));
            gridvec = cat(2, tpsGridX(:), tpsGridY(:))';
            
            % generate the surface
            rxX = rij(gridvec,X);
            GxX = rxX.^2.*log(rxX);
            GxX(rxX==0)=0;

            tpsGridZ = GxX*a + b(1) + b(2)*gridvec(1,:)' + b(3)*gridvec(2,:)';
            tpsGridZ = reshape(tpsGridZ, Ngrid(1), Ngrid(2));

            % assign result
            this.fittedParam.tpsGridX = tpsGridX;
            this.fittedParam.tpsGridY = tpsGridY;
            this.fittedParam.tpsGridZ = tpsGridZ;
        end
        
    end
end