classdef RegionOfInterest < handle_light
    % Define Cartesian coordinate system and coordinate ranges in 3d, this
    % is used in PointCloud to define a transformation between the
    % frame of detected points and the raw data
    %
    % ROI.translation, rotation take the basis (and origin) into the 
    % stack coordinate system:
    %
    % translation*[originp 1]' = [0 0 0 1]
    % rotation*[xp 1]' = [1 0 0 1]
    %
    % which means points are moved to the stack coordinates with the
    % inverse, e.g. translation^{-1} takes the origin [0,0,0] to the stack
    % coordinates of that origin, originp
    
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

        rotation        % affine matrix, rotation into the new coordinates
        translation     % affine matrix, translation into the new coordinates
        
        xpRange         % [ xpmin xpmax ]
        ypRange         % [ ypmin ypmax ]
        zpRange         % [ zpmin zpmax ]
    end
    
    properties (Dependent)
        
        originp % the new origin provided by the centroid
        
        xp  	% new x axis in stack frame
        yp      % new y axis in stack frame
        zp      % new z axis in stack frame

        theta   % theta, phi - polar coordinates of the zp axis
        phi     % theta, phi - polar coordinates of the zp axis
        
        % chi - azimuthal angle to rotate x and y by after having rotated 
        %       the whole system first through theta at phi=0, then through phi  
        chi
        irrep   % minimal data to specify ROI, used in experiment class
    end
    
    %---------------------------------------------------------------------
    % protected methods used in visualization
    %---------------------------------------------------------------------
    
    methods (Static = true) %, Access = protected
        
        function [i, j] = view2ind(view)
            % VIEW2IND Convert view name to indices for first and second axis
            %
            % [i, j] = view2ind(view)
            % 
            % e.g. view2ind('xy') returns [1 2]
            
            if strcmp(view,'xy')
                i = 1;  j = 2;
            elseif strcmp(view, 'yz')
                i = 2;  j = 3;
            elseif strcmp(view, 'zx')
                i = 3;  j = 1;
            end
        end
    end
    
    methods (Access = protected)
        
        function vert = getBoxVertices(this, xpRange, ypRange, zpRange)
            % GETBOXVERTICES Get vertex coordinates for a box specified by
            %  ranges in the primed coordinates
            %
            % vert = getBoxVertices(xpRange, ypRange, zpRange)
            %
            % vert is an 8x3 array whose rows are the corners of the box
            
            % vertices can be labeled as a 3 bit number: xyz are binary
            vert = zeros([8 3]);
            for k = 1:8
                ind = bitget(k-1, 3:-1:1) + 1;
                vert(k,:) = this.originp - [1 1 1]...
                            + xpRange(ind(1))*this.xp...
                            + ypRange(ind(2))*this.yp...
                            + zpRange(ind(3))*this.zp;
            end
        end
	end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    
    
    methods
    
        function this = RegionOfInterest(rotation, translation)
            % REGIONOFINTEREST Create a new RegionOfInterest
            %
            % RegionOfInterest(rotation, translation)
            %
            % rotation, translation should be affine (4x4) matrices
            
            this.rotation = rotation;
            this.translation = translation;
        end
        
        %--------------------------------------------------------
        % getters for dependent vars
        %--------------------------------------------------------
        
        function xp = get.xp(this)
            % get x coordinate in primed system
            xp = inv(this.rotation)*[1 0 0 0]';
            xp = xp(1:3)';
        end
        
        function yp = get.yp(this)
            % get y coordinate in primed system
            yp = inv(this.rotation)*[0 1 0 0]';
            yp = yp(1:3)';
        end
        
        function zp = get.zp(this)
            % get z coordinate in primed system
            zp = inv(this.rotation)*[0 0 1 0]';
            zp = zp(1:3)';
        end
        
        function originp = get.originp(this)
            % get origin coordinate in primed system
            originp = inv(this.translation)*[0 0 0 1]';
            originp = originp(1:3)';
        end
        
        function theta = get.theta(this)
            theta = acos(dot([0 0 1], this.zp)); 
        end
        
        function phi = get.phi(this)
            phi = atan2(this.zp(2), this.zp(1));
        end

        function chi = get.chi(this)
            % rotation that takes zp -> z
            R1 = rotationmat3D(-this.theta, [0 1 0])*rotationmat3D(-this.phi, [0 0 1]);
            xpp = R1*this.xp';
            chi = atan2(xpp(2), xpp(1));
        end
        
        function irrep = get.irrep(this)
            % the independent properties that should be saved by experiment
            
            irrep = struct( 'rotation', this.rotation,...
                            'translation', this.translation,...
                            'xpRange', this.xpRange,...
                            'ypRange', this.ypRange,...
                            'zpRange', this.zpRange);
        end
                
        %--------------------------------------------------------
        % setters
        %--------------------------------------------------------
        
        function setAngles(this, phi, theta, chi)
            % SETANGLES Convert angles into coordinate axes, set axes and origin
            %
            % setAngles(phi, theta, chi)
            %
            % theta, phi:   polar coordinates of the zp axis
            % chi:          azimuthal angle to rotate x and y by after 
            %               having rotated the whole system first through 
            %               theta at phi=0, then through phi
            
            rotationtmp = rotationmat3D(-phi, [0 0 1]);
            rotationtmp = rotationmat3D(-theta, [0 1 0])*rotationtmp;
            this.rotation(1:3,1:3) = rotationmat3D(-chi, [0 0 1])*rotationtmp;
        end
        
        function setOrigin(this, originp)
            % SETORIGIN Set the originp property
            % 
            % setOrigin(originp)
            %
            % originp:  the ROI origin relative to the stack
            
            % set originp    
            this.translation = [1 0 0 -originp(1);...
                                0 1 0 -originp(2);
                                0 0 1 -originp(3);
                                0 0 0 1];
            % also update ranges
            this.xpRange = this.xpRange - originp(1);
            this.ypRange = this.ypRange - originp(2);
            this.zpRange = this.zpRange - originp(3);
        end
        
        function setAxes(this, xp, yp, zp)
            % SETAXES Set the ROI coordinate axes
            %
            % setAxes(xp, yp, zp)
            %
            % xp, yp, zp should be row vectors
            % other variables are automatically updated to stay consistent
            
            % first find the angles of zp
            theta = acos(dot([0 0 1], zp)); 
            phi = atan2(zp(2), zp(1));
            
            % rotation1 rotates zp into z
            rotation1 = rotationmat3D(-theta, [0 1 0])*rotationmat3D(-phi, [0 0 1]);
            % and xp into xpp
            xpp = rotation1*xp';
            % chi is the angle required to rotate xpp into x
            chi = atan2(xpp(2), xpp(1));
            %rotationmat3D(-chi, [0 0 1])*xpp;
            
            % then set those
            this.setAngles(phi, theta, chi);
        end 
        
        function setRanges(this, xpRange, ypRange, zpRange)
            % SETRANGES Set ranges defining the actual region of interes
            %
            % setRanges(xpRange, ypRange, zpRange)
            %
            % xpRange, ypRange zpRange: ranges in the coordinate system 
            %                           defined by xp, yp, zp, originp
            
            if ~((numel(xpRange) == 2 && numel(ypRange) == 2 && numel(zpRange) ==2)...
              || (numel(xpRange) == 0 && numel(ypRange) == 0 && numel(zpRange) ==0))
                error('range is specified as [xmin xmax] or [] empty');
            end
            
            this.xpRange = xpRange;
            this.ypRange = ypRange;
            this.zpRange = zpRange;
        end
        
        function setAlignment(this, R, T)
            % SETALIGNMENT Set the alignment of the ROI
            %
            % setAlignment(R,T)
            % 
            % R, T: affine matrices specifying rotation and translation
            %
            % this is an alternative to setting the origin, axes or angles
            
            this.rotation = R;
            this.translation = T;
        end
        
        function newROI = composeAlignment(this, rotation, translation)
            % COMPOSEALIGNMENT Compose current alignment with new relative alignment
            %
            % composeAlignment(rotation, translation)
            %
            % rotation, translation: homogeneous matrices
            %
            % this does not update xpRange etc, needs to be done
            % externally!
            
            % let a be the translation vector, then R*T(a) = R + T(Ra) - 1
            % so rotation^(-1) acting on the last column of R*T will give a,
            % we can decompose R*T*R'*T' into a translation and rotation
            % the same way
            
            oldRT = rotation*translation;
            newRT = this.rotation*this.translation;
            totalRT = newRT*oldRT;
            
            totalR = eye(4);
            totalR(1:3,1:3) = totalRT(1:3,1:3);

            totalT = eye(4);
            totalT(:,4) = inv(totalR)*totalRT(:,4);

            newROI = surfaceDetection.RegionOfInterest(totalR, totalT);
        end

        % ------------------------------------------------------
        % visualization
        % ------------------------------------------------------
        
        function drawAxes(this, view, s)
            % DRAWAXES Draw coordinate axes projected on xy,yz or zx plane
            %
            % drawAxes(view, s)
            %
            % view: string specifying view, for example 'xy'
            % s:    scale, number specifying axis length
            
            O = this.originp;
            xp = this.xp; yp = this.yp; zp = this.zp;
            
            [i, j] = this.view2ind(view);
            
            lw = 2;     % line width
            
            line([O(i), O(i) + s*xp(i)], [O(j), O(j) + s*xp(j)],...
                                            'Color', 'r', 'LineWidth', lw);
            line([O(i), O(i) + s*yp(i)], [O(j), O(j) + s*yp(j)],...
                                            'Color', 'g', 'LineWidth', lw);
            line([O(i), O(i) + s*zp(i)], [O(j), O(j) + s*zp(j)],...
                                            'Color', 'b', 'LineWidth', lw);
        end
             
        function drawROI(this, view)
            % DRAWROI Draw ROI projected on the xy,yz or zx plane
            %
            % drawROI(view)
            %
            % view: string specifying view, for example 'xy'
            
            
            % vertices of the ROI
            vert = this.getBoxVertices(this.xpRange, this.ypRange, this.zpRange);

            % edges are lines between 4 special vertices and the vertices
            % that are one bit flip away from them
            svert = [bin2dec('000') bin2dec('011') bin2dec('110') bin2dec('101')] + 1;

            ec = 'm';   % edge color
            lw = 1;     % line width
            
            [i, j] = this.view2ind(view);
            
            for k = svert(1:end)
                for m = 1:3
                    % flip a bit
                    n = bitset(k-1, m, ~bitget(k-1,m)) + 1;
                    % draw an edge
                    line([vert(k, i) vert(n, i)], [vert(k, j) vert(n, j)],...
                        'Color', ec, 'LineWidth', lw);
                end
            end
        end

    end
end