classdef radialEdgeDetector < surfaceDetection.surfaceDetector
    % Surface detection based on thresholded radial derivative with pre-
    % and post processing to get a cleaner result
    %
    % First reduces the background, then computes cylindrical radial Gaus- 
    % sian derivative. Instead of finding the maximal intensity jump, 
    % a user provided threshold is used to generate a binary mask. 
    % This mask is cleaned in each plane along z by removing small groups 
    % of isolated pixels and a morphological closure and the (x,y) position
    % of pixels on the perimeter is added to the point cloud.
    %
    % options:  channel : channel of the stack to use
    %           sigma : width of the Gaussian 
    %           ssfactor : sub-sampling factor
    %           rmRadialOutliers : remove radial outlier >
    %               rmRadialOutliers*sigma (set zero for none)
    %           rmIntensityOutliers : remove edge intensity outliers
    %           thresh : threshold of intensity to be foreground
    %           amin : minimal connectec foreground pixel size
    %           bgdisc : radius of disc for background estimation
    %           dildisc : radius of dilation disc
    %           sp: stack dimension permutation     
    
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
            'rmRadialOutliers', 2, 'rmIntensityOutliers', 2,...
            'thresh',25,'amin',20,'bgdisc',0,'dildisc',20,'sp',[1 2 3]);      
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = radialEdgeDetector()
            % Constructor
            %
            % radialEdgeDetector()
            
            this = this@surfaceDetection.surfaceDetector();
        end 
        
        % ------------------------------------------------------
        % surface detection
        % ------------------------------------------------------
        
        function detectSurface(this, stack)
            % Detect surface in the stack with preset options.
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
            
            debugMsg(1, ['radialEdge.detectSurface() : channel='...
                num2str(opts.channel) ', sigma=' num2str(opts.sigma)...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', rmRadialOutliers=' num2str(opts.rmRadialOutliers)...
                ', rmIntensityOutliers=' num2str(opts.rmIntensityOutliers)...
                ', threshold =' num2str(opts.thresh)...
                ', min area =' num2str(opts.amin)...
                ', background = ' num2str(opts.bgdisc) ...
                ', dilation disc =' num2str(opts.dildisc) '\n']);
            
            if ~isempty(stack)
                data = stack.image.apply{opts.channel};
            else
                error('stack is empty');
            end
            
            if sum( (opts.sp-[1 2 3]).^2 )~=0
                % permute the stack; 
                data = permute(data,opts.sp([2 1 3]));
            end
            
            if opts.ssfactor > 1
                data = data(1:opts.ssfactor:end, 1:opts.ssfactor:end, 1:opts.ssfactor:end);
            else
                if numel(data) > 10^8
                    crazy = questdlg('Are you sure you want to call the detector with no sub-sampling?',...
                        'Big data, no sub-sampling, potentially exceeding available memory, slow progress!',...
                        'Yes','No','No');
                    if strcmp(crazy, 'No')
                        return;
                    end
                end
            end
            
            
            %---------------------------------
            % edge detection 
            %---------------------------------

            debugMsg(2, 'Edge detection: ');
            ticID = tic;

            %debugMsg(2, 'Applying Laplacian of Gaussian, ');
            %LOGdata = abs(LOG(data, opts.sigma));
            
            % background Correction
            if opts.bgdisc > 0
                
                se = strel3D('sphere', opts.bgdisc);
                dataEroded = imerode(data, se);
                dataErodedDilated = imdilate(dataEroded,se);
                data = data-dataErodedDilated;
            end
            % edge detection
            debugMsg(2, 'cylindrical radial derivative');
            edgeData = -cylinderRadialD(data, opts.sigma, 1);
            dt = toc(ticID);
            debugMsg(2, ['dt = ' num2str(dt) ' sec\n']);
            
            %---------------------------------
            % point cloud detection
            %---------------------------------
            
            %---------------------------------
            % zero pad edge data
            %---------------------------------
            ps = 4;
            padsize = opts.dildisc+ps;

            ypadding = zeros([padsize, size(edgeData,2), size(edgeData, 3)], class(edgeData));
            edgeData = cat(1, ypadding, edgeData, ypadding);

            xpadding = zeros([size(edgeData,1), padsize, size(edgeData, 3)], class(edgeData));
            edgeData = cat(2, xpadding, edgeData, xpadding);

            zpadding = zeros([size(edgeData,1), size(edgeData, 2), padsize], class(edgeData));
            edgeData = cat(3,zpadding, edgeData, zpadding);
            
            if sum( (opts.sp-[1 2 3]).^2 )~=0
                
                inverseorder(opts.sp) = 1:numel(opts.sp);
                ySize = size(edgeData,inverseorder(2));
                xSize = size(edgeData,inverseorder(1));
                zSize = size(edgeData,inverseorder(3));
            else 
                ySize = size(edgeData,1);
                xSize = size(edgeData,2);
                zSize = size(edgeData,3);
            end
            seD = strel('disk',opts.dildisc);
            
            zmin = 1;
            zmax = size(edgeData, 3);
            
            bwseg = edgeData < -1000;
            
            pointCloud = zeros(0,3);
            for z = zmin:zmax

                d = edgeData(:,:,z)>opts.thresh;
                d = bwareaopen(d,opts.amin);
                d = imdilate(d,seD);
                d = imfill(d,'holes');
                d(1,:) = 0;
                d(end,:) = 0;
                d(:,1) = 0;
                d(:,end) = 0;
                d = bwareaopen(imerode(d,seD),opts.amin);
                bwseg(:,:,z) = d;
            end
            
            % separated determination of pointcloud from above, to single
            % out small connected components.
            bwseg = bwlabeln(bwseg);
            si = regionprops(bwseg);
            Areas = cat(1,si.Area);
            [~,ii] = max(Areas);
            if ~isempty(ii)
            bwseg(bwseg~=ii(1)) = 0;
            for z = zmin:zmax
               
                per = bwperim(bwseg(:,:,z));
                [I,J] = ind2sub(size(per),find(per));
                ZZ = z*ones(size(I));
                % don't forget to remove padding offset from above;
                pointCloud = [pointCloud;J-opts.dildisc-ps,I-opts.dildisc-ps,ZZ-opts.dildisc-ps];
            end
            else
                pointCloud = zeros(1,3);
            end
            %---------------------------------
            % remove radial outliers
            %---------------------------------
            redPointCloud = [];
            if opts.rmRadialOutliers > 0
                
                for z = zmin:zmax
                    
                    pcSlice = pointCloud(pointCloud(:,3) == z,:);
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
            
            %---------------------------------
            % scale point cloud to full size
            %---------------------------------
            
            largePC = zeros(size(pointCloud));
            for i = 1:3
                largePC(:,i) = (pointCloud(:,i)-1)*opts.ssfactor + 1;
            end
            
            largePC = sortrows(largePC, 3);
            if sum( (opts.sp-[1 2 3]).^2 )~=0
                % permute back; 
                inverseorder(opts.sp) = 1:numel(opts.sp);
                largePC = largePC(:,inverseorder([2 1 3]));
            end
            
            ROI = surfaceDetection.RegionOfInterest(eye(4),eye(4));
            % ranges
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
