function SOImask = SOImask3D(mesh, onionOpts, maskSize)
    % this function was used to check that 2D methods correctly sum
    % intensity 
    %
    % SOImask = SOImask3D(mesh, onionOpts, maskSize)
    %
    % SOImask:  binary mask
    % 
    % mesh:         structure with field f, v, vn for faces, vertices,
    %               vertex normals
    % onionOpts:    onion options as provided to
    %               SurfaceOfInterest.pullbackStack
    % maskSize:     required dimensions of mask
    %
    % requires VOXELIZE from matlab file exchange

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
    
    dN = onionOpts.layerDistance*(onionOpts.nLayers-1)/2;

    normal = mesh.vn./repmat(sqrt(sum(mesh.vn.^2,2)),[1 3]);
    Vouter = mesh.v + dN*normal;
    Vinner = mesh.v - dN*normal;

    meshFV = struct();
    meshFV.faces = mesh.f;

    meshFV.vertices = Vouter;
    maxes = round(max(meshFV.vertices));
    mins = round(min(meshFV.vertices));
    sizes = maxes - mins;
    SOImaskOuter = VOXELISE(sizes(1),sizes(2),sizes(3),meshFV);

    meshFV.vertices = Vinner;
    maxes = round(max(meshFV.vertices));
    mins = round(min(meshFV.vertices));
    sizes = maxes - mins;
    SOImaskInner = VOXELISE(sizes(1),sizes(2),sizes(3),meshFV);

    offset = round(min(Vinner)) - round(min(Vouter));
    bla = SOImaskOuter(offset(1):offset(1)+sizes(1)-1, offset(2):offset(2)+sizes(2)-1,...
                    offset(3):offset(3)+sizes(3)-1);
    preSOImask = SOImaskOuter;
    preSOImask(offset(1):offset(1)+sizes(1)-1, offset(2):offset(2)+sizes(2)-1,...
                    offset(3):offset(3)+sizes(3)-1) = bla - SOImaskInner;

    SOImask = false(maskSize);
    maxes = round(max(Vouter));
    mins = round(min(Vouter));
    sizes = maxes - mins;
    SOImask(mins(2):mins(2)+sizes(2)-1, mins(1):mins(1)+sizes(1)-1,mins(3):mins(3)+sizes(3)-1) = permute(preSOImask, [2 1 3]);
end