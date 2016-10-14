classdef PointCloud < handle_light
    % Contains the point cloud that represents a detected surface
    %
    % property ROI specifies the coordinate system and range in which the
    % points are specified relative to the the raw data 
    %
    % See also surfaceDetection.RegionOfInterest
       
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
        points;             % Nx3 matrix containing point coordinates
        ROI;                % RegionOfInterest object
    end
    
    properties (Dependent)
        unalignedPoints;    % point coordinates in raw data reference frame
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    
    methods
    
        function this = PointCloud(points, varargin)
            % POINTCLOUD Create a new PointCloud
            %
            % PointCloud(points)
            % PointCloud(points, ROI)
            %
            % points:   Nx3 array containing point coordinates
            % ROI:      RegionOfInterest object specifying coordinates
            %
            % When passing ROI, it doesn't transform the points, i.e.
            % assumes that the points are in the aligned frame, to
            % transform, call constructor without alignment and use setROI
            
            this.points = points;
            
            % alignment
            if nargin == 2
                this.ROI = varargin{1};
            else
                this.ROI = surfaceDetection.RegionOfInterest(eye(4), eye(4));
            end
            
            % default ranges
            this.determineRanges(0);
        end
        
        % ------------------------------------------------------
        % point cloud moments
        % ------------------------------------------------------
        
        function [centroid, cov] = getMoments(this)
            % GETMOMENTS Get centroid and central second moments of point cloud
            %
            % [centroid, cov] = getMoments()

            if isempty(this.points)
                error('pointCloud is empty');
            end
            
            D = 3;
            centroid = mean(this.points,1);
            centered = zeros(size(this.points));
            for i=1:D
                centered(:,i) = double(this.points(:,i)) - centroid(i);
            end
            
            cov = zeros(D);
            for i = 1:D
                for j = 1:D
                    summand = centered(:,i).*centered(:,j);
                    cov(i,j) = mean(summand(:));
                end
            end
        end
        
        % ------------------------------------------------------
        % point cloud alignment
        % ------------------------------------------------------
        
        function relativeAlign = determineROI(this, margin)
            % DETERMINEROI Set ROI with from point cloud
            %
            % relativeROI = determinROI(margin)
            %
            % margin: margin of bounding box around the point cloud
            %
            % alignment based on the principle axes and the centroid 
            % boundingbox (xpRange, ..) on points in this frame
            %
            % ROI has the alignment properties:
            %   - xp yp zp: the new coordinate axes in the old frame these are 
            %   the principal axes unless the eigenvalues are near degenerate
            %   - originp : the new origin provided by the centroid
            %   - R, T: rotation and translation the new coordinates
            %   - theta, phi, psi : the angles that define xp, yp, zp
            %
            % relativeAlign is the ROI being set relative to the current
            % this is useful when updating multiple pointclouds
            %
            % See also surfaceDetection.PointCloud.determineRanges
            
            % calculate moments for alignment
            [originp, M] = this.getMoments();

            x = [1 0 0]';
            z = [0 0 1]';

            % major axis is z-axis, the others are arbitrary
            [V, D]= eig(M);
            D = diag(D);

            % sort eigenvalues by size
            [~,idx] = sort(D, 'descend');

            % if largest ev not much larger than smallest, do nothing
            if ( D(idx(1)) - D(idx(3)) )/( D(idx(1)) + D(idx(3)) ) < 0.1
                disp('all components are similar, returning identity');
                % don't return yet. In this case we need to also update the
                % roi!
                
                relativeAlign = surfaceDetection.RegionOfInterest(eye(4), eye(4));
                relativeAlign.setAxes(this.ROI.xp,this.ROI.yp,this.ROI.zp);
                relativeAlign.setOrigin(originp);
            
                % the current alignment should be composed with the relative
                % alignment found here
            
                this.setROI(relativeAlign);
            
                this.determineRanges(margin);
                return;
            % else make the principle axis the z-axis
            else
                zp = V(:, idx(1));
            end
            zp = real(zp);
            % if the remaining eigenvalues are not very different, keep x and y in the
            % same planes
            if ( D(idx(2)) - D(idx(3)) )/( D(idx(2)) + D(idx(3)) ) < 0.1
                debugMsg(2, 'smaller eigenvalues are similar, keeping x and y in plane\n');
                yp = cross(zp, x);
                yp = yp/norm(yp);
                xp = cross(yp, zp);
            % else make the second axis thee y axis
            else
                yp = V(:, idx(2));
                xp = cross(yp, zp);
            end
            
            relativeAlign = surfaceDetection.RegionOfInterest(eye(4), eye(4));
            relativeAlign.setAxes(xp',yp',zp');
            relativeAlign.setOrigin(originp);
            
            % the current alignment should be composed with the relative
            % alignment found here
            
            this.setROI(relativeAlign);
            this.determineRanges(margin);
        end
        
        function determineRanges(this, margin)
            % DETERMINERANGES Boundingbox (xpRange, ..) based on point cloud
            %
            % determineRanges(margin)
            %
            % margin: margin of bounding box around the point cloud
            %
            % See also surfaceDetection.PointCloud.determineROI
            
            maxs = ceil(max(this.points, [], 1));
            mins = floor(min(this.points, [], 1));

            xpRange = [mins(1)-margin, maxs(1)+ margin];
            ypRange = [mins(2)-margin, maxs(2)+ margin];
            zpRange = [mins(3)-margin, maxs(3)+ margin];
            
            this.ROI.setRanges(xpRange, ypRange, zpRange);
        end
            
        %------------------------------------------------------
        % get points in stack frame
        %------------------------------------------------------
        
        function unalignedPoints = get.unalignedPoints(this)
            % Get the pointcloud in the stack frame

            Rinv = inv(this.ROI.rotation);
            Tinv = inv(this.ROI.translation);
            
            alignedPoints = [this.points ones([size(this.points,1) 1])];
            % (RT)^(-1) = T^-1 R^-1
            unalignedPoints = (Tinv*Rinv*alignedPoints')';
            unalignedPoints = unalignedPoints(:,1:3);
        end
        
        % ------------------------------------------------------
        % set alignment and ROI
        % ------------------------------------------------------
         
        function setROI(this, relativeROI)
            % SETROI Set alignment (composes), transform points accordingly
            %
            % setROI(relativeROI)
            %
            % relativeROI:  RegionOfInterest object specifying axes
            %               relative to this ROI
            
            R = relativeROI.rotation;
            T = relativeROI.translation;
            
            pointsHomogenous = [this.points ones([size(this.points,1) 1])];
            alignedPoints = (R*T*pointsHomogenous')';
            alignedPoints = alignedPoints(:,1:3);
            
            this.points = alignedPoints;

            % compose relative roi with old roi ( = this.ROI)            
            this.ROI = relativeROI.composeAlignment(this.ROI.rotation, this.ROI.translation);
            this.determineRanges(10);
        end
        
        % ------------------------------------------------------
        % set points
        % ------------------------------------------------------
        
        function setPoints(this, points)
            % SETPOINTS Set the point coordinate matrix
            %
            % setPoints(points)
            %
            % points:   3xN matrix with point coordinates
            
            this.points = points;
        end
        
        % ------------------------------------------------------
        % visualize point cloud
        % ------------------------------------------------------
        
        function inspect(this, ssfactor)
            % INSPECT 3D plot of point cloud
            %
            % inspect(ssfactor)
            %
            % ssfactor: subsampling factor of the point cloud
            
            ind = 1:ssfactor:length(this.points);
            dotsize = 4;
            scatter3(this.points(ind,1), this.points(ind,2), this.points(ind,3),...
                dotsize, this.points(ind,3), 'filled');
            axis equal; 
            colormap Jet
        end
        
    end
end