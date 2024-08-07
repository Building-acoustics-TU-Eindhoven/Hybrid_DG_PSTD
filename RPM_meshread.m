function [VX,VY,Nv,K,EToV,BClines,Ekind] = RPM_meshread(msh_file)

mesh = load_gmsh2(msh_file);

if ~exist(msh_file,'file')
    error('mesh file does not exist!')
end

VX1 = mesh.POS(:,1);
VY1 = mesh.POS(:,2);
Nv = length(VX1);

if mesh.nbTriangles & mesh.nbQuads
    error('Elements of mixed types are not possible (TRANGLES+QUADS')
end

if mesh.nbTriangles
    K = mesh.nbTriangles;
    EToV_1 = mesh.TRIANGLES([1:K],[1:3]);
    Ekind = 'triangles';
end
if mesh.nbQuads
    K = mesh.nbQuads;
    EToV_1 = mesh.QUADS([1:K],[1:4]);
    Ekind = 'quads';
end

BClines1 = mesh.LINES([1:mesh.nbLines],:);

% Fix mesh
[~,EToV,ridx]=RPM_fixmesh([VX1 VY1],EToV_1);
VX = VX1(ridx).';
VY = VY1(ridx).';
% modify boundary lines accordingly
BClines = zeros(size(BClines1));
for ii=1:size(BClines1,1)
    id1 = find(ridx==BClines1(ii,1));
    id2 = find(ridx==BClines1(ii,2));
    BClines(ii,:) = [id1 id2 BClines1(ii,3)];
end
Nv = length(VX);
K = length(EToV(:,1));
Nv = length(VX);
Nel_DG = K;

% Reorder elements to ensure orientation
if strcmp(Ekind,'quads')
    ax = VX(EToV(:,1));
    ay = VY(EToV(:,1));
    bx = VX(EToV(:,2));
    by = VY(EToV(:,2));
    cx = VX(EToV(:,3));
    cy = VY(EToV(:,3));
    dx = VX(EToV(:,4));
    dy = VY(EToV(:,4));

    D_quads = (ax.*by + bx.*cy + cx.*dy + dx.*ay - bx.*ay - cx.*by - dx.*cy - ax.*dy);
    indx_quads = find(D_quads<0);
    EToV(indx_quads,:) = EToV(indx_quads,[1 4 3 2]);
elseif strcmp(Ekind,'triangles')
    % not implemented
    warning('orientation of geometrical nodes for triangles might be wrong')
end


end
