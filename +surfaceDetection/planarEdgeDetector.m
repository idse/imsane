classdef planarEdgeDetector < surfaceDetection.planarDetector
    % planarDetector Extract a surface point cloud from a planar data stack
    %
    % As opposed to other detectors, the detected surface is
    % represented not only by a PointCloud object but also by an image
    % surfaceMatrix, containing z values for each xy
    
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
        
        maxdIdx %maximal derivative in z
    end
    
    properties (Constant)
        % default detector options
        defaultOptions = struct('sigma', 1, 'channels', 1, 'zdir', 3,...
                                'maxIthresh', 0, 'summedIthresh', 0,...
                                'sigZoutliers', 1, 'scaleZoutliers', 3,...
                                'seedDistance', 20); 
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = planarEdgeDetector()
            % PLANARDETECTOR Create a new planarDetector
            %
            % planarDetector()
            
            this = this@surfaceDetection.planarDetector();
        end 
        
        % ------------------------------------------------------
        % surface detection
        % ------------------------------------------------------
        
        function detectSurface(this, stack, seed)
            % DETECTSURFACE detects the surface as the position of the 
            % maximal Gaussian derivative in some direction z, i.e. the position of the
            % largest intensity jump along some direction and smoothened over some
            % scale
            %
            % options:  sigma :     width of the Gaussian z-derivative
            %           channels :  channels (summed) to use for detection
            %           zdir :      dimension corresponding to z, minus flips direction
            %
            %           maxIthresh:     throw out points with MIP dimmer than this
            %           summedIthresh:  throw out points with SIP dimmer than this
            %
            %           both maxIthresh and summedIthresh can be specified
            %           as scaled values between 0 and 1 or absolute values
            %           between 1 and 2^n-1 where where n is hte number of
            %           bits
            %
            %           sigZoutliers:   remove height outliers after all other masks
            %           scaleZoutliers: spatial scale of outlier removal
            %
            % options are set using setDetectOptions
            % the options in the second block only affect filters / masks on the
            % already detected surface and can be modified without redetecting
            %
            % scaleZoutliers is the linear size of a region over which the
            % distribution of height is computed, sigZoutliers is then a cutoff in
            % units of standard deviation of this distribution to remove misdetected
            % points far above or below the other points in the region
            %
            % See also surfaceDetection.surfaceDetector.setOptions
           
            % default seed is empty
            if nargin == 2
                seed = [];
            end
            
            opts = this.options;
            
            if ~isfield(opts, 'channels')
                disp('detectSurface: options.channel not specified, using channel 1');
                opts.channels = 1;
            end
            
            if ~isfield(opts, 'zdir')
                opts.zdir = 3;
            end
            
            % now add up all the indicated channels for detection
            imageGrids = stack.image.apply;
            tdata = zeros(size(imageGrids{1}), class(imageGrids{1}));
            
            w = stack.Ilim{1}(2)./[stack.Ilim{:}];
            w = w(2:2:end);
            
            for i = 1:length(opts.channels)
                tdata = tdata + w(opts.channels(i)).*imageGrids{opts.channels(i)};
            end
            
            zSize = size(tdata, 3); 
            xSize = size(tdata, 2);
            ySize = size(tdata, 1);
            
            %---------------------------------
            % stack mask
            %---------------------------------

            % determine stack mask boundingbox to save time on the
            % detection in the next step
            
            mask = stack.getMask();
            
            % for some region regionprops boundingbox manages to make my
            % computer run out of memory, so  
            
            % xz mask
            projMask = squeeze(max(mask,[],1))';
            [X,Z] = meshgrid(1:xSize, 1:zSize);
            xmin = min(X(projMask));
            xmax = max(X(projMask));
            zmin = min(Z(projMask));
            zmax = max(Z(projMask));

            % yz mask
            projMask  = squeeze(max(mask,[],2))';
            [Y,Z] = meshgrid(1:ySize, 1:zSize);
            ymin = min(Y(projMask));
            ymax = max(Y(projMask));
            zmin = max(min(Z(projMask)), zmin);
            zmax = min(max(Z(projMask)), zmax);

            % xymask
            projMask = squeeze(max(mask,[],3));
            [X,Y] = meshgrid(1:xSize, 1:ySize);
            xmin = max(min(X(projMask)), xmin);
            xmax = min(max(X(projMask)), xmax);
            ymin = max(min(Y(projMask)), ymin);
            ymax = min(max(Y(projMask)), ymax);
            
            xmin = max(floor(xmin), 1);
            xmax = min(ceil(xmax), xSize);
            ymin = max(floor(ymin), 1);
            ymax = min(ceil(ymax), ySize);
            zmin = max(floor(zmin), 1);
            zmax = min(ceil(zmax), zSize);
            
            %---------------------------------
            % intensties
            %---------------------------------
            
            this.summedI = sum(tdata(ymin:ymax, xmin:xmax, zmin:zmax), 3);
            this.maxI = max(tdata(ymin:ymax, xmin:xmax, zmin:zmax), [], 3);
            
            %---------------------------------
            % edge detection 
            %---------------------------------
            
            debugMsg(2, 'Edge detection: ');
            ticID = tic;
            
            % convolve with a Gaussian kernel and insert zeros where the 
            % convolution can't be obtained without zero padding the original image
            
            if size(5*opts.sigma) > zSize
                error(['Gaussian kernel size ' num2str(size(kernel))...
                    ' greater than image size ' num2str(zSize)]);
            end
            
            dIdx = tdata(ymin:ymax, xmin:xmax, zmin:zmax);
            
            for j = 1:3
                % first derivative if j == 3 (z direction)
                % otherwise zeroth derivative: smoothing
                dIdx = myconvn(dIdx, -single(GaussD(opts.sigma, j==3, j)));
            end

            % just take positive part and rescale for easy display
            if sign(opts.zdir) > 0
                dIdx(dIdx < 0) = 0;
            else
                dIdx(dIdx > 0) = 0;
                dIdx = -dIdx;
            end

            % make a logical array of the maxima of the derivative
            maxima = max(dIdx, [],3);
            maxima3d = repmat(maxima, [1, 1, zmax - zmin + 1]);
            apical = (dIdx == maxima3d & dIdx > 0.01);

            this.maxdIdx = maxima;
            
            % to turn this into a mesh, take the maximal z coordinate for fixed x,y
            % slice by slice to avoid using too much memory
            [~,Z] = meshgrid(xmin:xmax, zmin:zmax);
            
            this.surfaceMatrix = zeros([ySize xSize]);

            for y = ymin:ymax
                this.surfaceMatrix(y,xmin:xmax) = ...
                                    max(squeeze(apical(y-ymin+1,:,:)).*Z', [], 2);
            end
            
            % apply the Z-projection of the stack mask
            surfMatMask = max(mask,[],3);
            this.surfaceMatrix = surfMatMask.*this.surfaceMatrix;
            
            dt = toc(ticID);
            debugMsg(2, ['dt = ' num2str(dt) ' sec\n']);
            
            %---------------------------------
            % point cloud detection
            %---------------------------------
            
            debugMsg(2, 'Detecting point cloud, ');
            ticID = tic;
            
            % generate the point cloud
            Z = this.surfaceMatrix;
            [X,Y] = meshgrid(1:xSize, 1:ySize);
            points = [X, Y, Z];
            this.pointCloud = surfaceDetection.PointCloud(points);
            this.pointCloud.ROI.setRanges([1 xSize], [1 ySize], [1 zSize]);
            
            % apply masks to clean up the point cloud
            this.applyMasks();
            
            % continuity based on seed
            %---------------------------
            if ~isempty(seed)
                
                dist = abs(this.surfaceMatrix - seed);
                this.mask = this.mask & (dist < opts.seedDistance);
                
                this.updatePointCloud();
            end
            
            dt = toc(ticID);
            debugMsg(2, ['dt = ' num2str(dt) ' sec\n']);
        end
    end
end