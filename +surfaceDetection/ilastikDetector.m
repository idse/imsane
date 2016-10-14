classdef ilastikDetector < surfaceDetection.surfaceDetector
    % Segmentation based on prediction maps from ilastik. 
    % 
    % 
    % options:  channel : channel of the stack to use
    %           sigma : width of the Gaussian 
    %           ssfactor : sub-sampling factor used in the ilastik
    %           classifier.
    %           rmRadialOutliers : remove radial outlier >
    %               rmRadialOutliers*sigma (set zero for none)
    %           thresh : threshold of intensity to be foreground between 0
    %                                   and 1
    %           amin : minimal connectec foreground pixel size
    %           dildisc : radius of dilation disc  
    %           fileName: Name of the ilastik prediction file. Typically
    %               ending with Probabilities.h5 in Ilastik v1.1
    %           foreGroundChannel: the labeled used in ilastik as
    %               foreground default is 2. (for two classes, 1 would then be
    %               the background).
    %           zdim: cylinder axis in matlab coords, 2 = x
    
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
        defaultOptions = struct('channel', 1, 'sigma', 2, 'ssfactor', 4,...
            'rmRadialOutliers', 2,'thresh',.5,'amin',20,'dildisc',10,...
            'fileName',[],'foreGroundChannel',2,'zdim',2);      
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = ilastikDetector()
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
            
            debugMsg(1, ['ilastik.detectSurface() : channel='...
                num2str(opts.channel) ', sigma=' num2str(opts.sigma)...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', rmRadialOutliers=' num2str(opts.rmRadialOutliers)...
                ', threshold =' num2str(opts.thresh)...
                ', min area =' num2str(opts.amin)...
                ', dilation disc =' num2str(opts.dildisc)...
                ', fileName =' num2str(opts.fileName)...
                ', fileGroundChannel =' num2str(opts.foreGroundChannel)...
                ', zDim =' num2str(opts.zdim),'\n']);
            
            if isempty(opts.fileName)
                error('Please provide a regular prediction from ilastik in h5 format.');
            end

            %---------------------------------
            % Segmentation of a prediction map from ilastik. 
            %---------------------------------
            
            
            % load the exported data out of the ilastik prediction
            fileName = [opts.fileName,'_Probabilities.h5'];
            h5fileInfo = h5info(fileName);
            if strcmp(h5fileInfo.Datasets.Name,'exported_data')
                file = h5read(fileName,'/exported_data');
            elseif strcmp(h5fileInfo.Datasets.Name,'volume')
                file = h5read(fileName,'/volume/prediction');
            else
                error(['Please provide a regular prediction from ilastik, either in', ...
                    'the format of version 1.1 or 0.5']);
            end
            
            foreGround = opts.foreGroundChannel;
            
            % ilastik internally swaps axes. 1: class, 2: y, 3: x 4 : z 
            pred = permute(file,[3,2,4,1]);
            pred = pred(:,:,:,foreGround);
            pred = uint8(255*pred);
           
            % size of prediction
            idxPerm = circshift(1:3, [1 -opts.zdim]);

            ySize = size(pred, idxPerm(1));
            xSize = size(pred, idxPerm(2));
            zSize = size(pred, idxPerm(3));
            
            zmin = 1;
            zmax = size(pred, idxPerm(3));
            
            % loop through each slize along the 2nd dimension and perform
            % thresholding with morphological operations to obtain a closed
            % curve.
            points = struct('x',[],'y',[],'z',[]);
            for z = 1 : zSize

                
                if opts.zdim == 3
                    im = squeeze(pred(:,:,z));
                elseif opts.zdim == 2
                    % transpose needed to preserve chirality!
                    im = squeeze(pred(:,z,:));
                else
                    im = squeeze(pred(z,:,:));
                end

                mask = bwareaopen(im > 255*opts.thresh,opts.amin);
                
                mask([1,end],:) = 0;
                mask(:,[1,end]) = 0;
                mask = imdilate(mask,strel('disk',opts.dildisc));
                mask = imfill(mask,'holes');
                seg  = imerode(mask,strel('disk',opts.dildisc));

                per = bwperim(seg);
                if z > 1
                ind = find(per==1);
                else
                    ind = find(seg == 1);
                end
                [I,J] = ind2sub(size(seg),ind);

                % extract the point cloud. 
                points(z).x = I;
                points(z).y = J;
                points(z).z = ones(size(I))*z;

            end

            % don't forget to rescale the data;
            x = cat(1,points.x);
            y = cat(1,points.y);
            z = cat(1,points.z);     
            
            
            pointCloud = [x,y,z];
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
            if opts.zdim == 2 % i.e zp == x bc x is the second index
                ROI.setAxes([0 1 0], [0 0 1], [1 0 0]);
            elseif opts.zdim == 1    
                ROI.setAxes([0 0 1], [1 0 0], [0 1 0]);
            end
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
        
        % ------------------------------------------------------
        % prepare data for ilastik segmentation
        % ------------------------------------------------------
        
        function prepareIlastik(this,stack)
            % Prepaere stack for Ilastik segmetnation.
            %
            % prepareIlastik(stack)
            
            % Accoring to the specified options sub-sample the stack and
            % save it for analysis with ilastik. 
            
            opts = this.options;
            
            im = stack.image.apply();
            
            fileName = [opts.fileName,'.h5'];
            
            if exist(fileName,'file');
                delete(fileName)
            end    
                  
            dsetName = '/inputData';
            
            
            for c = 1 : length(im)
 
                %for k = 1 : opts.ssfactor: size(im{c},3)
                    
                %    image = im{c}(:,:,k);
                %    imwrite(image(1:opts.ssfactor:end,1:opts.ssfactor:end),...
                %        fileName,'tiff','Compression','none','WriteMode','append');
                %end
                
                image(:,:,:,c) = im{c}(1:opts.ssfactor:end,1:opts.ssfactor:end,1:opts.ssfactor:end);
            end
            if ndims(image)==4
                h5create(fileName,dsetName,[size(image,2) size(image,1), size(image,3) size(image,4)]);
                h5write(fileName,dsetName,permute(image,[2 1 3 4]));
            else
                h5create(fileName,dsetName,[size(image,2) size(image,1), size(image,3)]);
                h5write(fileName,dsetName,permute(image,[2 1 3]));
            end
            
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
