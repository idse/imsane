classdef pointCloudDetector < surfaceDetection.surfaceDetector
    % This is not a detector, it rather reads a previously detected point 
    % cloud and assembles it into a detector object. 
    % 
    % The point cloud should be Npoints rows and 3 columns. CSV or txt are
    % currently supported. 
    % 
    % options: 
    %           ssfactor : sub-sampling factor used when the point cloud
    %                       was detected. Assumes isotropic downsampling. 
    %           rmRadialOutliers : remove radial outlier >
    %               rmRadialOutliers*sigma (set zero for none)
    %           fileName: Name of the ilastik prediction file. Typically
    %               ending with Probabilities.h5 in Ilastik v1.1
    %           zdim: cylinder axis in matlab coords, 2 = x
    %           readerType: csv or dlm reader.
    %           yInvert: Matlab inverts y-axes, by default this option is
    %               1. turn it off if the detection seems swapped.
    
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2014 Idse Heemskerk and Sebastian Streichan
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
    
    properties (Constant)
        % default detector options
        defaultOptions = struct('ssfactor', 4,'rmRadialOutliers', 2,...
            'fileName',[],'zdim',2,'readerType','csv','yInvert',1);      
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = pointCloudDetector()
            % Constructor
            %
            % radialEdgeDetector()
            
            this = this@surfaceDetection.surfaceDetector();
        end 
        
        % ------------------------------------------------------
        % surface detection
        % ------------------------------------------------------
        
        function detectSurface(this, stack)
            % Detect surface will read the point cloud from disc and
            % assmble a detector object.
            %
            % detectSurface(stack)
            
            % Detect the surface using cylindricalRadialD function from
            % generalFunctions/filters with specified options.
            % First down-sample the stack, then corect for background, 
            % detect edges, slice by slice compare intensity to thresh and
            % apply filters specified in the options. 
            % With resulting point cloud correct for radial outliers
            % and intensity outliers out. 
            
            opts = this.options;
            
            debugMsg(1, ['pointCloudDetector.detectSurface() :'...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', rmRadialOutliers=' num2str(opts.rmRadialOutliers)...         
                ', fileName =' num2str(opts.fileName)...
                ', zDim =' num2str(opts.zdim)...
                ', readerType =' opts.readerType...
                ', yInvert = ' num2str(opts.yInvert),'\n']);
                
            
            if isempty(opts.fileName)
                error(['Please provide a regular file name for the point cloud.'...
                    'The format should be CSV (comma separated values) or ASCI' ... 
                    'delmited file supported by matlab. See doc dlmread or doc csvread.']);
            end

            
            %---------------------------------
            % load a csv file from disc
            %---------------------------------
            
            
            % load the exported data out of the ilastik prediction
            fileName = opts.fileName;
            
            if strcmp(opts.readerType,'csv')
                pointCloud = csvread(fileName);
            elseif strcmp(opts.readerType,'dlm')
                pointCloud = dlmread(fileName,'\t');
                %pointCloud = pointCloud(:,[2,1,3]);
            else
                error(['Please provide a regular file name for the point cloud.'...
                    'The format should be CSV (comma separated values) or ASCI' ... 
                    'delmited file supported by matlab. See doc dlmread or doc csvread.']);
            end
            
            %---------------------------------
            % remove radial outliers
            %---------------------------------
            
            %idxPerm = circshift(1:3, [1 -opts.zdim]);

            ySize = stack.imageSize(1);
            xSize = stack.imageSize(2);
            zSize = stack.imageSize(3);
            
            if opts.yInvert == 1
                pointCloud(:,2) = xSize-pointCloud(:,2);
            end
            
            zmin = 1;
            zmax = zSize;
            
            
            %pointCloud = pointCloud(:,idxPerm);
            
            redPointCloud = [];
            if opts.rmRadialOutliers > 0
                
                for z = zmin:zmax
                    
                    pcSlice = pointCloud(pointCloud(:,opts.zdim) == z,:);
                    CM = mean(pcSlice, 1);
                    Rsq = zeros([size(pcSlice,1) 1]);
                    for i = 1:3
                        Rsq = Rsq + (double(pcSlice(:,i)) - CM(i)).^2;
                    end
                    R = sqrt(Rsq);
                    good =  abs(R - mean(R)) < opts.rmRadialOutliers*std(R);
                    pcSlice = pcSlice(good, :);
                    redPointCloud = [redPointCloud;pcSlice];
                end
            end
            pointCloud = redPointCloud;
            
            
            %--------------------------------------------------
            % scale point cloud to full size and set alignment
            %--------------------------------------------------
            
            largePC = zeros(size(pointCloud));
            for i = 1:3
                largePC(:,i) = (pointCloud(:,i)-1)*opts.ssfactor + 1;
            end
            
            largePC = sortrows(largePC, 3);
            
            debugMsg(2, 'setting ROI alignment based on zDim\n');
            
            ROI = surfaceDetection.RegionOfInterest(eye(4), eye(4));
            %if opts.zdim == 2 % i.e zp == x bc x is the second index
            %    ROI.setAxes([0 1 0], [0 0 1], [1 0 0]);
            %elseif opts.zdim == 1    
            %    ROI.setAxes([0 0 1], [1 0 0], [0 1 0]);
            %end
            
            
            % we also want to set the ranges of the ROI, i.e. some bounding 
            % box that contains the pointcloud which will be used by the
            % fitter to determine its initial domain
            xpRange = [1, xSize*opts.ssfactor];
            ypRange = [1, ySize*opts.ssfactor];
            zpRange = [1, zSize*opts.ssfactor];
            ROI.setRanges(xpRange,ypRange,zpRange);
            this.pointCloud = surfaceDetection.PointCloud(largePC,ROI);
            this.pointCloud.determineROI(15);
            %this.pointCloud = surfaceDetection.PointCloud(largePC);

        end
        
        
        function inspectQuality(this, inspectOpts, stack)
            %   inspect quality of fit in single slice in dimension specified 
            %   by options and display image.
            %
            %   inspectQuality(inspectOpts, stack)
            %
            %   override inspectQuality to deal with subsampled point cloud
            
            ssfac = this.options.ssfactor;
            
            if ssfac ~= 1 && rem(inspectOpts.value-1, ssfac) ~= 0
                inspectOpts.value = round(inspectOpts.value/ssfac)*ssfac + 1;
                disp([]);
                debugMsg(2, ['WARNING: sub-sampled point cloud taken from different'...
                    ' plane, may look a little off\n']);
            end

            inspectQuality@surfaceDetection.surfaceDetector(this, inspectOpts, stack);
        end
        
    end
end
