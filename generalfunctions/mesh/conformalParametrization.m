function uv = conformalParametrization(mesh, constrVidx, constrV)
    % compute conformal parametrization
    %
    % SCP minimizes a functional u'*L_C*u + lambda*(u'*B*u-1)
    % where the first term is the conformal energy with L_C the Laplacian
    % minus the the area matrix 
    % and the second adds a constraint u'*B*u where B is the boundary
    % matrix
    %
    % we generalize this by replace B by K, an arbitrary kernel that in
    % particular can enforce temporal continuity in a sequence of meshes
    %
    % based on the implementation of SCP by Ryan Schmidt

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
    
    if nargin == 1
        constrVidx = [];
        constrV = [];
    end
    
    % faces and vertices
    faces = mesh.f;
    verts = mesh.v;
    Nverts = size(verts,1);

    % % list of boundary vertex indices
    bidx = mesh.b{1};
    
    % vertex face one ring
    vfring = compute_vertex_face_ring(faces);

    % construct initial NxN weight matrix
    W = -compute_mesh_weight(verts, faces, 'conformal');
    W = -W; % negate so that sums on diagonal are positive
    W(1:Nverts+1:Nverts*Nverts) = -sum(W,2);  % set W(i,i) = -sum_j W(i,j)

    % need to construct 2Nx2N system that includes both X and Y 
    % uv-values. We will use form [X,0; 0,Y]
    Ld = [W,sparse(Nverts,Nverts); sparse(Nverts,Nverts), W];

    % area matrix
    A = sparse([],[],[],2*Nverts,2*Nverts,4*numel(bidx));
    for ii = 1:numel(bidx)
        i = bidx(ii);    ix = i;    iy = ix+Nverts;

        vtris = vfring{i};

        for ti = 1:numel(vtris)
            f = faces(vtris(ti), :);
            [j,k] = tripick(f,i);
            jx = j;  jy = j+Nverts;
            kx = k;  ky = k+Nverts;

            A(ix,ky) = A(ix,ky) - 1;
            A(ix,jy) = A(ix,jy) + 1;

            A(iy,jx) = A(iy,jx) - 1;
            A(iy,kx) = A(iy,kx) + 1;
        end
    end
    
    % construct natural conformal system L_C   [Mullen08]
    Lc = Ld-A;
    
    % if no constraints are provided, 
    % default to spectral conformal parametrization
    if isempty(constrV)
        
        Bi = [bidx, bidx+Nverts];
        B = sparse(Bi,Bi,ones(numel(Bi),1),2*Nverts,2*Nverts);
        
        % solve eigensystem
        eigsoptions.disp = 0; 
        eigsoptions.isreal = 1; 
        eigsoptions.issym = sum(sum(abs(Lc-Lc'))) < eps;
        neigs = 5;
        fudge = 10e-8;      % numerical fudge from sec 3.4
        [V,D] = eigs(Lc + fudge*speye(2*Nverts,2*Nverts), B, neigs, 0, eigsoptions);
        %    [V,D] = eigs(B, Lc + fudge*speye(2*N,2*N), neigs, 'lm', eigsoptions);
        k = 3;
        uv = real([V(1:Nverts, k), V(Nverts+1:2*Nverts, k)]);  
        
        
    % else add the constraints to the matrix equation and solve 
    else
        
        constrU = constrV(:,1);
        constrV = constrV(:,2);

        Nc = 2*numel(constrVidx);
        
        constraints = zeros([Nc 2*Nverts]);

        for i = 1:Nc/2
            constraints(i, constrVidx(i)) = 1;
            constraints(Nc/2 + i, Nverts + constrVidx(i)) = 1;
        end

        up = zeros([2*Nverts + Nc, 1]);
        up(2*Nverts + (1:Nc/2)) = constrU;
        up(2*Nverts + (Nc/2+1:Nc)) = constrV;

        Lcp = cat(1, Lc, constraints);
        V = mldivide(Lcp, up);

        uv = real([V(1:Nverts), V(Nverts+1:2*Nverts)]);  
        
    end
end