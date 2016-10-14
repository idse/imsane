classdef (Abstract) surfaceFitter < handle_light
    %   Fit a surface to a point cloud 
    %
    %   This is an abstract class that needs to be implemented for
    %   specific cases. For examples see the tpsFitter or spherelikeFitter
       
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

        fitDomain   % FiniteSet describing the domain on which the fit is done
        
        fittedPoints  % cell array {X,Y,Z} of points on the fitted surface
                      % depending on fitter XYZ may have to be grids
                      
        % this replaces in a way compatible with meshes: 
        % embedding   % embedding patch needed for normally evolve, alignment, inspect quality
        
        % for example: order of polynomal fit, defined per subclass
        %
        % every fitter has these options:
        %
        % shift:        shift by which to normally evolve after fitting
        % subsampling:  surface subsampling in normally evolve
        fitOptions  = struct(...        
            'shift', [],...             
            'normEvolveSS', []);        
   
        fittedParam % parameters as result of fitting
    end
    
    properties (Abstract, SetAccess = protected)
       
        % charts - charts this fitter can produce and their properties
        %
        % the rows have the structure: 
        % {name, description, stepsize, intersections, fundamental, desired}
        %
        % stepsize:         2d vector of stepsize
        % intersections:    lists indices of charts with overlapping domain
        % fundamental:      boolean, does chart define a set in the topology
        % desired:          boolean indicating whether to produce it
        charts               
    end

    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------

    methods (Abstract)
        
        % FITSURFACE Fit the pointcloud to produce an embedding patch
        %
        % Implementation should set fitDomain and fittedParam.
        % Fit options can be set through through the generic function
        % setFitOptions. The subclass defines what can be set.
        %
        % see also surfaceFitting.surfaceFitter.setFitOptions
        fitSurface(this);
       
        % GENERATEEMBEDDING Create the embedding in some chart
        %
        % generateEmbedding(chartName)
        %
        % The embedding is generated using the result of a fit, stored
        % in fittedParam, so fitSurface should be called first.
        % The available charts are listed in the charts property of the
        % SurfaceFitter object.
        % 
        % See also surfaceFitting.surfaceFitter.fitSurface
        generateEmbedding(this, chartName);    
       
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
        % See also surfaceFitting.surfaceFitter.fitSurface
        populateSOI(this, SOI);                 
        
        % NORMALLYEVOLVE Normally shift surface out or in by some amount
        %
        % normallyEvolve(shift)
        %
        % shift : pixel distance along the normal
        normallyEvolve(this, shift);         
    end
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = surfaceFitter()
            % SURFACEFITTER Create a new SurfaceFitter
            %
            % surfaceFitter()
            %
            % Initializes default fitOptions
            % Should also create empty fittedParam structure in subclasses

            this.fitOptions.shift = 0;
        end
        
        % ------------------------------------------------------
        % setters
        % ------------------------------------------------------
        
        function setChartResolution(this, chartName, stepSize)
            % SETCHARTRESOLUTION Set the step size: resolution of charts
            %
            % setChartResolution(chartName, stepSize)
            %
            % chartName:     name of a chart
            % stepSize:      vector with 2 entries

            % determine the number of charts
            ncharts = length(this.charts);
            % loop through charts and set stepsize for the matching
            % chart
            for k = 1 : ncharts
                if strcmp(chartName, this.charts(k).name)
                    this.charts(k).stepSize = stepSize;
                end
            end
        end
    
        function setFitOptions(this, fitOptions)
            % SETFITOPTIONS Set options used in the specific fit
            %
            % setFitOptions(fitOptions)
            %
            % fitOptions:   structure containing the options, available
            % options are defined in each subclass of surfaceFitter, see
            % e.g. spherelikeFitter
            %
            % See also surfaceFitting.spherelikeFitter.fitSurface
            
            fields = fieldnames(this.fitOptions);
            
            for i = 1:numel(fields)
                if isfield(fitOptions, fields{i})
                    this.fitOptions.(fields{i}) = fitOptions.(fields{i});
                else
                    debugMsg(1, ['fitOptions.' fields{i} ' not provided, stays the same\n']);
                end
            end
        end
       
        function setFittedParameters(this, fittedParam)
            % SETFITTEDPARAMETERS Set parameters obtained from the fit
            %
            % fittedParam:  structure with entries corresponding to the
            %               fits obtained depending on the fit used. 
            %
            % These are in principle determined by fitSurface, but
            % sometimes it is useful to set them by hand, e.g. when
            % reloading a previously made fit.
            %
            % See also surfaceFitting.spherelikeFitter
            
            fields = fieldnames(this.fittedParam);
            
            for i = 1:numel(fields)
                if isfield(fittedParam, fields{i})
                    this.fittedParam.(fields{i}) = fittedParam.(fields{i});
                else
                    debugMsg(1, ['fittedParam.' fields{i} ' not provided, stays the same\n']);
                end
            end
        end
        
        function setFitDomain(this, fitDomain)
            % SETFITDOMAIN Set domain on which to fit. 
            %
            % setFitDomain(fitDomain)
            %
            % fitDomain:    finiteSet specifying the boundary, stepsize
            %               gridsize, alignedEmbeddingSpace and the dimension.
            %
            % See also diffgeometry.FiniteSet
            
            this.fitDomain = fitDomain;
        end
        
        % ------------------------------------------------------
        % setter for desired charts
        % ------------------------------------------------------
        
        function setDesiredChart(this, chartName, desired)
            % SETDESIREDCHART Set which charts should be computed. 
            %
            % setDesiredChart(chartName, desired)
            %
            % chartName:    string matching one of the available charts
            % desired:      boolean, 0 for no and 1 for yes.
            
            debugMsg(1, 'surfaceFitter.setDesiredChart()\n');
            
            found = false;
            
            for i = 1 : length(this.charts)
                
                availChartName = this.charts(i).name;
                if strcmp(chartName, availChartName)
                    found = true;
                    this.charts(i).desired = desired;
                end
            end
            
            if ~found
                debugMsg(1, ['chart ' chartName ' does not exist\n']);
            end
        end
        
        
        % ------------------------------------------------------
        % getters
        % ------------------------------------------------------

        function stepSize = getChartResolution(this, chartName)
            % GETCHARTRESOLUTION Get stepsize of fundamental chart
            %
            % getChartResolution(chartName)
            %
            % chartName:    string corresponding to the name of the chart
            
            ncharts = length(this.charts);
            for k = 1 : ncharts
                if strcmp(chartName,this.charts(k).name)
                    stepSize = this.charts(k).stepSize;
                end
            end
           
        end
        
        % ------------------------------------------------------
        % Write Quality inspection to disc.
        % ------------------------------------------------------
        
        function storeQualityInspection(this, options, varargin)
            % STOREQUALITYINSPECTION Inspect fit quality and store to file
            %
            % storeQualityInspection(options, detector, stack)
            %
            % options:              struct with fields:
            % - inspectOptions:     options for inspectQuality
            % - range:              range of planes along dimension
            % - outName:            filename of stored figure
            % - export:             boolean write to file
            % - closeFig:           boolean close figure after
            %
            % See also surfaceFitting.surfaceFitter.inspectQuality
            
            if ~isfield(options,'range')   
                options.range = 1:this.domain.boundary{3}(2);
            end
            if ~isfield(options,'outName');
                options.outName = 'qualityPointCloudFit.tif';
            end
            if ~isfield(options,'export');
                options.export = 'true'; 
            end
            if ~isfield(options,'closeFig')
                options.closeFig = 'true';
            end
            if ~isfield(options,'inspectOptions')
                options.inspectOptions = struct();
            end
            
            if strcmp(options.closeFig,'true')
                h = figure;
            end
          
            % loop through all optionally set values (default is full size)
            for value = options.range;
                options.inspectOptions.value = value;
                this.inspectQuality(options.inspectOptions, varargin{:})
                
                if strcmp(options.export,'true')
                    if (value-options.range(1))>0
                        I = getframe(h);
                        imwrite(I.cdata, options.outName,'tiff','Compression','none','WriteMode','append');
                    else
                        I = getframe(h);
                        imwrite(I.cdata, options.outName,'tiff','Compression','none');
                    end
                end
            end

            if strcmp(options.closeFig,'true')
                close(h);
            end
           
        end
        
        % ------------------------------------------------------
        % inspect quality of single slice and display image.
        % ------------------------------------------------------
        
        function inspectQuality(this, options, pointCloud, stack)
            % INSPECTQUALITY Inspect quality of fit in specified cross section
            %
            % inspectQuality(options, detector, stack)
            % inspectQuality(options, pointcloud, stack)
            %   
            % options:          structure with the following fields:
            %   - dimension:    'x', 'y' or 'z', normal to the cross section
            %   - value:        position of cross section along dimension
            %   - pointCloud:   color in which to display point cloud, e.g.
            %                   'b', to not display, remove this option
            %   - noalign       dont rotate fitted points according to
            %                   detector coordinate system
            % pointCloud        PointCloud or (for compatibility)
            %                   SurfaceDetector
            % stack:            Stack object (containing raw data)
           
            % pointCloud 
            if isa(pointCloud, 'surfaceDetection.surfaceDetector')
                pointCloud = pointCloud.pointCloud;
            end
            
            % parse the display options
            if ~isfield(options,'dimension'), options.dimension = 'z'; end
            if ~isfield(options,'noalign'), options.noalign = false; end
            
            if options.dimension == 'x' 
                ind = [1 2 3];
            elseif options.dimension == 'y'
                ind = [2 3 1];
            else
                ind = [3 1 2];
            end

            X = this.fittedPoints;
            alignedFitPts = [X{1}(:),X{2}(:),X{3}(:)];
            
            % show image, possibly with pointcloud
            imshow(stack.getSlice(options.dimension, options.value)); 
            
            if ~isempty(pointCloud)
                
                lspec = ['*' options.pointCloud];
                PC = pointCloud.unalignedPoints;
                
                PCplane = PC(round(PC(:,ind(1))) == options.value, :);
                % ugly way to get some points from the upsampled cloud
                if isempty(PCplane)
                    PCplane = PC(abs(PC(:,ind(1)) - options.value) < 1, :);
                end
                hold on 
                plot(PCplane(:,ind(2)), PCplane(:,ind(3)), lspec, 'MarkerSize', 2)    
                hold off
                
                % take detector alignment/ROI
                fitPC = surfaceDetection.PointCloud(alignedFitPts, pointCloud.ROI);
            else
                
                disp(['surfaceFitter.inspectQuality: empty detector passed,'...
                    ' not displaying detected point cloud']);
                % no alignment from detector
                % is ther an alignment from the xp? 
                fitPC = surfaceDetection.PointCloud(alignedFitPts);
            end
            
            % plot fit over the image 
            if options.noalign
                fitPoints = fitPC.points;
            else
                fitPoints = fitPC.unalignedPoints;
            end
            fitPlane  = fitPoints(round(fitPoints(:,ind(1))) == options.value, :);
            
            % a hack to add more points if the fitted points are sparse
            if numel(fitPlane) < 100
                d = 5;
                fitPlane  = fitPoints(round(fitPoints(:,ind(1))) <= options.value + d ...
                                        & round(fitPoints(:,ind(1))) >= options.value - d, :);
            end
            
            hold on 
            plot(fitPlane(:,ind(2)), fitPlane(:,ind(3)), 'ow', 'LineWidth', 2, 'MarkerSize', 2, 'Color', 'r')    
            hold off
        end
    end
end
