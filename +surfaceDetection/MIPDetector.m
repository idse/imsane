classdef MIPDetector < surfaceDetection.planarDetector
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
        
        function this = MIPDetector()
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
            % maximal intensity in some direction z
            %
            % seed:     height matrix (e.g. Z-component of embedding(t-1))
            %
            % options:  sigma :     width of the Gaussian smoothing before
            %           channels :  channels (summed) to use for detection
            %           zdir :      dimension corresponding to z, minus flips direction
            %
            %           maxIthresh:     throw out points with MIP dimmer than this
            %           summedIthresh:  throw out points with SIP dimmer than this
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
            for i = 1:length(opts.channels)
                tdata = tdata + imageGrids{opts.channels(i)};
            end
            
            zSize = size(tdata, 3); 
            xSize = size(tdata, 2);
            ySize = size(tdata, 1);
            
            %---------------------------------
            % point cloud detection
            %---------------------------------
            
            debugMsg(2, 'Detecting point cloud, ');
            ticID = tic;
            
            % smoothing
            if size(3*opts.sigma) > zSize
                error(['Gaussian kernel size ' num2str(size(kernel))...
                    ' greater than image size ' num2str(zSize)]);
            end            
            
            for j = 1:3
                % tdata = myconvn(tdata, single(GaussD(opts.sigma, 0, j)));
                tdata = convn(tdata, single(GaussD(opts.sigma, 0, j)), 'same');
            end

            % find positions of maximal intensity
            this.summedI = sum(tdata, 3);
            [this.maxI, this.surfaceMatrix] = max(tdata, [], abs(opts.zdir));
            
            % adjust maxI
            this.maxI = imadjust(mat2gray(this.maxI));

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