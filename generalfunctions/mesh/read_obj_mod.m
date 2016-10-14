function mesh = read_obj_mod(filename)
% read_obj - load a .obj file.
%
%   mesh = read_obj(filename);
%
%   mesh.v :    vertices
%   mesh.f :    faces
%   mesh.vn:    vertex normals
%
%   Copyright (c) 2008 Gabriel Peyre
%
% modified by Idse Heemskerk 2015, to accomodate obj format as output by
% meshlab, return mesh struct, and include vertex normals

fid = fopen(filename);
if fid<0
    error(['Cannot open ' filename '.']);
end

frewind(fid);
vertex = [];
faces = [];
normal = [];
while 1
    s = fgetl(fid);
    if ~ischar(s), 
        break;
    end
    if ~isempty(s) && strcmp(s(1), 'f')
        % face
        %faces(:,end+1) = sscanf(s(3:end), '%d %d %d');
        bla = sscanf(s(3:end), '%d//%d %d//%d %d//%d');
        faces(:,end+1) = bla(1:2:end);
    end
    if numel(s) > 1 && strcmp(s(1:2), 'v ')
        % meshlab separates decimals by commas on some systems
        s = strrep(s, ',', '.'); 
        % vertex
        vertex(:,end+1) = sscanf(s(3:end), '%f %f %f');
    end
    if numel(s) > 1 && strcmp(s(1:2), 'vn')
        % vertex normal
        s = strrep(s, ',', '.');
        normal(:,end+1) = sscanf(s(3:end), '%f %f %f');
    end
end
faces = faces';
vertex = vertex';
normal = normal';

mesh = struct('v', vertex, 'f', faces, 'vn', normal);

fclose(fid);

