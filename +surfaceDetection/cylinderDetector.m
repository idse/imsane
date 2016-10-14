classdef cylinderDetector < surfaceDetection.surfaceDetector
    % Detect surface by maximal cylindrical radial derivative
    % as opposed to fast cylinder detector, the center for the radial
    % derivative is estimated separately for each z-position based on
    % thresholding
   
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
    
    properties (Constant)
        % default detector options
        defaultOptions = struct('channel', 1, 'sigma', 2, 'ssfactor', 4,...
            'nBins', 120, 'rmRadialOutliers', 2, 'rmIntensityOutliers', 2,...
            'zDim', 3);     
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = cylinderDetector()
            % FASTCYLINDERDETECTOR Create a new FastCylinderDetector
            %
            % fastCylinderDetector()
            
            this = this@surfaceDetection.surfaceDetector();
        end 
        
        % ------------------------------------------------------
        % surface detection
        % ------------------------------------------------------
        
        function detectSurface(this, stack)
            % DETECTSURFACE Detect the surface using cylindrical radial
            % derivative
            %
            % detectSurface(stack)
            %
            % options: channel:             channel of the stack to use
            %          sigma:               width of the Gaussian 
            %          ssfactor:            sub-sampling factor
            %          nBins:               number of angular bins
            %          rmRadialOutliers:    remove radial outlier > rmRadialOutliers*sigma 
            %                               (set to zero for none)
            %          rmIntensityOutliers: remove edge intensity outliers
            %          zDim :               cylinder axis in matlab coords: 2 = x
            %
            % options are set using setDetectOptions
            %
            % See also surfaceDetection.surfaceDetector.setOptions
            
            
            % First down-sample the stack, then detect edges, slice by
            % slice extract resulting point cloud taking radial outliers
            % and intensity outliers out. 
            opts = this.options;
            
            debugMsg(1, ['fastCylinder.detectSurface() : '...
                'channel=' num2str(opts.channel)...
                ', sigma=' num2str(opts.sigma)...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', zDim=' num2str(opts.zDim)...
                ', rmRadialOutliers=' num2str(opts.rmRadialOutliers)...
                ', rmIntensityOutliers=' num2str(opts.rmIntensityOutliers) '\n']);
            
            if ~isempty(stack)
                data = stack.image.apply{opts.channel};
            else
                error('stack is empty');
            end
            
            % downsample for detection, 
            % (don't convert to 8bit, there is no speed gain and little
            % memory)
            if opts.ssfactor > 1
                data = data(1:opts.ssfactor:end, 1:opts.ssfactor:end, 1:opts.ssfactor:end);
            else
                if numel(data) > 10^8
                    crazy = questdlg('Are you sure you want to call the detector with no sub-sampling?',...
                        'Big data, no sub-sampling, potentially a bad idea',...
                        'Yes','No','No');
                    if strcmp(crazy, 'No')
                        return;
                    end
                end
            end
%             
%             %---------------------------------
%             % edge detection 
%             %---------------------------------
%             
%             debugMsg(2, 'Edge detection: ');
%             ticID = tic;
%             
%             debugMsg(2, 'cylindrical radial derivative');
%             edgeData = -cylinderRadialD(data, opts.sigma, 1, opts.zDim);
%             
%             dt = toc(ticID);
%             debugMsg(2, ['dt = ' num2str(dt) ' sec\n']);
%             
            %---------------------------------
            % point cloud detection
            %---------------------------------
            
            debugMsg(2, 'Detecting point cloud, \n');
            
            ticID = tic;
            
            % third element of idxPerm is z dimension, first two are the remaining
            % dimensions
            idxPerm = circshift(1:3, [1 -opts.zDim]);
            
            ySize = size(data, idxPerm(1));
            xSize = size(data, idxPerm(2));
            zSize = size(data, idxPerm(3));
            yCenter = round(ySize/2); 
            xCenter = round(xSize/2); 

            % now get the entire pointcloud, slice by slice

            pointCloud = [];
            PCI = [];
            
            for z = 1:zSize
                
                % progress indicator
                fprintf('.');
                if rem(z,80) == 0
                    debugMsg(1,'\n');
                end
                
                % point cloud indices in slice
                pcSlice = zeros([opts.nBins 3], 'uint16');
                % correspond edge intensities (sharpness)
                PCISlice = zeros([opts.nBins 1]);
                
                if opts.zDim == 3
                    slice = squeeze(data(:,:,z));
                elseif opts.zDim == 2
                    % transpose needed to preserve chirality!
                    slice = squeeze(data(:,z,:))';
                else
                    slice = squeeze(data(z,:,:));
                end
                
                %center of object
                dslice = mat2gray(slice); % double needed for threshholding
                bwslice = dslice > graythresh(dslice);
                [y x] = find(bwslice);
                xCenter = round(mean(x));
                yCenter = round(mean(y));
                
                % radial edge detection
                edgeSlice = -cylinderRadialD(slice, opts.sigma, 1, opts.zDim, [xCenter yCenter]);
                
                % meshgrids need later
                [X, Y] = meshgrid(-xCenter:(xSize-xCenter)-1, -yCenter:(ySize-yCenter)-1);
                [theta] = cart2pol(X, Y);   
                                
                % set up the angle masks
                dphi = 2*pi/opts.nBins;
                anglemask = false([numel(theta) opts.nBins]);
                xs = {};
                ys = {};

                for n = 1:opts.nBins
                    x = -pi + n*dphi;
                    angleMaskGrid = (theta >= x) & (theta < x + dphi);

                    anglemask(:,n) = angleMaskGrid(:);
                    xs{n} = X(anglemask(:,n));
                    ys{n} = Y(anglemask(:,n));
                end

                %------

                for n = 1:opts.nBins
                    
                    % make a pixel mask for the maximal laplacian in a
                    % radial bin
                    edgeAngle = edgeSlice(anglemask(:,n));
                    maxEdgeI = min(edgeAngle);                  %HERE IS A CHANGE
                    maxEdgeIdx = (edgeAngle == maxEdgeI);
                    
                    % then find the coordinates for that maximal value
                    % if multiple pixels have the maximal value, average
                    if sum(maxEdgeIdx) == 1
                        pcSlice(n, 1) = xs{n}(maxEdgeIdx) + xCenter + 1;
                        pcSlice(n, 2) = ys{n}(maxEdgeIdx) + yCenter + 1;
                        PCISlice(n) = maxEdgeI;
                    else
                        pcSlice(n, 1) = mean(xs{n}(maxEdgeIdx)) + xCenter + 1;
                        pcSlice(n, 2) = mean(ys{n}(maxEdgeIdx)) + yCenter + 1;
                        PCISlice(n) = maxEdgeI;
                    end
                    pcSlice(n, 3) = z;
                end
                
                % throw out outliers in radius distribution of slice
                if opts.rmRadialOutliers > 0
                    
                    CM = mean(pcSlice, 1);
                    Rsq = zeros([size(pcSlice,1) 1]);
                    for i = 1:3
                        Rsq = Rsq + (double(pcSlice(:,i)) - CM(i)).^2;
                    end
                    R = sqrt(Rsq);
                    good =  abs(R - mean(R)) < opts.rmRadialOutliers*std(R);
                    pcSlice = pcSlice(good, :);
                    PCISlice = PCISlice(good, :);
                end
                
                % append points in slice to point cloud
                pointCloud = cat(1, pointCloud, pcSlice);
                PCI = cat(1, PCI, PCISlice);
            end
            fprintf('\n');

            dt = toc(ticID);
            debugMsg(2, ['dt = ' num2str(dt) ' secs\n']);
	
            % remove intensity outlier
            nSigI = opts.rmIntensityOutliers;
            if nSigI > 0
                pointCloud = pointCloud(PCI > mean(PCI) - nSigI*std(PCI), :);
            end
            
            %--------------------------------------------------
            % scale point cloud to full size and set alignment
            %--------------------------------------------------
            
            largePC = zeros(size(pointCloud));
            for i = 1:3
                largePC(:,i) = (pointCloud(:,i)-1)*opts.ssfactor + 1;
            end
            
            largePC = sortrows(largePC, 3);
            
            % the pointCloud comes out in the permuted coordinate system so
            % we have to pass this as an alignment to the constructor
            debugMsg(2, 'setting ROI alignment based on zDim\n');
            
            ROI = surfaceDetection.RegionOfInterest(eye(4), eye(4));
            if opts.zDim == 2 % i.e zp == x bc x is the second index
                ROI.setAxes([0 1 0], [0 0 1], [1 0 0]);
            elseif opts.zDim == 1    
                ROI.setAxes([0 0 1], [1 0 0], [0 1 0]);
            end
            
            % we also want to set the ranges of the ROI, i.e. some bounding 
            % box that contains the pointcloud which will be used by the
            % fitter to determine its initial domain
            xpRange = [1, xSize*opts.ssfactor];
            ypRange = [1, ySize*opts.ssfactor];
            zpRange = [1, zSize*opts.ssfactor];
            ROI.setRanges(xpRange, ypRange, zpRange);
            
            this.pointCloud = surfaceDetection.PointCloud(largePC, ROI);
        end
        
        function inspectQuality(this, inspectOpts, stack)
            % INSPECtQUALTIY Inspect quality of detection by overlaying a point cloud
            % through a slice in the data
            % 
            % inspectQuality(options, stack)
            %
            % inspect quality of single slice in dimension specified by options 
            % and display image.
            % options:  channels: array of channels onto RGB e.g. [1 2 NaN] 
            %           dimension: x,y,z 
            %           value
            %           pointCloud: color (matlab convention, g, w etc)
            %           ssfactor: subsampling of point cloud
            %
            % the last option is added in this subclass
            
            ssfac = this.options.ssfactor;
            
            if ssfac ~= 1 && rem(inspectOpts.value-1, ssfac) ~= 0
                inspectOpts.value = round(inspectOpts.value/ssfac)*ssfac + 1;
                debugMsg(2, ['WARNING: sub-sampled point cloud taken from different'...
                    ' plane, may look a little off\n']);
            end

            inspectQuality@surfaceDetection.surfaceDetector(this, inspectOpts, stack);
        end
        
    end
end
