function [subm, vidx, submFidx] = submeshFromV(mesh, vidx, clipears)
    % get a submesh by specifying vertex indices
    %
    % subm = submeshFromV(mesh, vidx)
    %
    % mesh : struct with fields v, f, vn, (optionally u)
    % subm : same
    % vidx : logical index

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
    
    if size(vidx,2) > 1
        vidx = vidx';
    end

    if nargin == 2
        clipears = true;
    end
    
    assert(islogical(vidx), 'vidx should be logical');
    
    % determine the corresponding faces in the submesh
    fidx = ismember(mesh.f, find(vidx));
    submFidx = all(fidx, 2);
    
    % it is possible that one of the points is not in any face, remove
    % those points
    newvidx = false(size(vidx));
    newvidx(unique(mesh.f(submFidx,:))) = true;
    vidx = newvidx;

%     % CLIPEARS PROBLEM WITH ATLAS BOUNDARIES< SOME INDEXING GOES WRONG
%     if clipears
%         
%         debugMsg(2, 'submeshFromV: clipping ears\n');
%         newf = mesh.f(submFidx,:);
%         [newf, newFidx] = clip_mesh(newf);
% 
%         submFidx(submFidx) = newFidx;
%         vidx = false(size(vidx));
%         vidx(unique(newf)) = true;
%     end
    
    % now create a proper submesh
    newv = mesh.v(vidx,:);
    newvn = mesh.vn(vidx,:);
    
    if isfield(mesh, 'u')
        newu = {};
        for i = 1:length(mesh.u)
            newu{i} = mesh.u{i}(vidx,:);
        end
    else
        newu = [];
    end

    old2new = zeros(size(vidx));
    old2new(vidx) = 1:sum(vidx);

    newf = mesh.f(submFidx,:);
    for j = 1:3
        newf(:,j) = old2new(newf(:,j));
    end

    % submesh boundaries
    b = compute_boundaries(newf);

    subm = struct('v', newv, 'f', newf, 'vn', newvn, 'u', {newu}, 'b', {b});
end
            