function strainStruct = metricStrainFromG(M)
    % get a strain structure from the metric at a point
    %
    % metricStrainFromG(M)
    %
    % M:    2x2 matrix
    
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
    
    
    linearstrain = eye(2) - (M + M')/2;
    rotation = (M(1,2) - M(2,1))/2;
    % V is the matrix of eigenvectors s.t. strain*V = V*D
    [V,D] = eig(linearstrain);

    % from this get the principle axis lengths and orientation in the same
    % format as region props
    % first figure out which eigenvalue is the larger one
    [n,~]=find(abs(D)==max(abs(D(:))));
    % the nondegenerate (generic case)
    if size(n,1)==1
        majoraxislength = D(n,n);
        minoraxislength = sum(D(D~=D(n,n)));
    %the degenerate case
    else
        n=1;
        majoraxislength = D(1,1);
        minoraxislength = D(2,2);
    end

    % if x of major axis is negative then rotate 180 degrees, so
    % orientation as in regionprops: angle from -90 (y) to 90 (y) degrees with x axis
    if V(1,n)<0
        major = -V(:,n);
    else
        major = V(:,n);
    end
    orientation = acos(major(1));
    % angle is negative if above the x axis
    if major(2) > 0
        orientation = -orientation;
    end
    orientation = orientation*180/pi;

    strainStruct = struct(  'MinorAxisLength', minoraxislength,...
                            'MajorAxisLength', majoraxislength,...
                            'Rotation', rotation,...
                            'Orientation', orientation );

end

