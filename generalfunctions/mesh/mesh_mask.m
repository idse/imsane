function mask = mesh_mask(mesh, uidx, minv, minu, vsize, usize)
% create a binary mask on a texture coordinate grid
%
% mask=mesh_mask(mesh, minv, minu, vsize, usize)
%
% mesh:         struct elements f,v,vn,u,b
% minv, minu:   offset of u,v grid
% vsize, usize: texture coordinate grid size

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

if ~isfield(mesh, 'b') || isempty(mesh.b)
    b = compute_boundaries(mesh.f);
else
    b = mesh.b;
end

assert(isfield(mesh, 'u') && ~isempty(mesh.u), 'u missing');

u = mesh.u{uidx};

% mask out interpolation outside the mesh 
nBdries = numel(b);

% determine which is the outer boundary
outBdry = 1;
for bi = 2:nBdries
    if all(inpolygon(u(b{1},1), u(b{1},2),...
                     u(b{bi},1),u(b{bi},2)));
        outBdry = bi;
    end
end

ub = u(b{outBdry}(:),1);
vb = u(b{outBdry}(:),2);

mask = poly2mask(ub - minu, vb - minv, vsize, usize);

% now cut out the holes
for bi = setdiff(1:nBdries, outBdry)

    ub = u(b{bi}(:),1);
    vb = u(b{bi}(:),2);

    inmask = poly2mask(ub -minu, vb - minv, vsize, usize);
    mask = mask & ~inmask;
end