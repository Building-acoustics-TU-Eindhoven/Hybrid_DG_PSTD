% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
Nfp = Norder+1; Np = (Norder+1)*(Norder+1); Nfaces=4;

% Compute nodal set
[r,s] = Nodes2D_quad(Norder);
r = r(:); s = s(:);

% Build reference element matrices
Vanderm = Vandermonde2D_quad(Norder,r,s); invVanderm = inv(Vanderm);
MassMatrix = invVanderm'*invVanderm;
[Dr,Ds] = Dmatrices2D_quad(Norder, r, s, Vanderm);

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)'; vd = EToV(:,4)';
x = zeros(length(r),Kel); y = zeros(length(r),Kel);
for k = 1 : Kel
    xc = VX(EToV(k,:)); yc = VY(EToV(k,:));
    [x(:,k),y(:,k)] = QuadMap(xc,yc,r(:),s(:));
end

% find all the nodes that lie on each edge
fmask1 = find( abs(s+1) < NODETOL)'; % south
fmask2 = find( abs(r-1) < NODETOL)'; % east
fmask3 = find( abs(s-1) < NODETOL)'; % north
fmask4 = find( abs(r+1) < NODETOL)'; % west
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

% Create surface integral terms
LIFT = Lift2D_quad(r,s,Fmask);

% calculate geometric factors
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);

% calculate geometric factors
[nx, ny, sJ] = Normals2D_quad(x,y,Dr,Ds,Fmask);
Fscale = sJ./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = tiConnect2D(EToV);

%% Allow for open boundary surfaces

%for ielem = 1:Kel
% TODO: get only elements potentially on boundary

% For each boundary surface element/face
for ibc = 1:size(BClines,1)

    % Get boundary vertices (2 vertices in 2D, 3 in 3D with tets)
    bcnode1 = BClines(ibc,1);
    bcnode2 = BClines(ibc,2);

    % Find elements sharing same vertices (ie boundary corresponds to element face)
    elem = [];
    for ielem = 1:Kel
        %for jj = 1:Nfaces
        if any(EToV(ielem,:)==bcnode1) && any(EToV(ielem,:)==bcnode2)
            % at least 1 element face corresponds to boundary face
            elem = [elem ielem];
        end
    end

    if length(elem)==2 % open surface, elem(1) and elem(2) on both sides of boundary

        % What if parallel mesh partition? Will we always see elements on
        % both side on open boundary?

        % Find common face between the two elements (corresponding to open boundary)
        face1 = find(EToE(elem(1),:)==elem(2)); % face # of element 1
        face2 = find(EToE(elem(2),:)==elem(1)); % face # of element 2

        % Check that face 'face1' of element 1 connects to face 'face2' of element 2
        assert(EToF(elem(1),face1)==face2)
        assert(EToF(elem(2),face2)==face1)

        % Manually remove connectivity
        EToE(elem(1),face1) = elem(1);
        EToF(elem(1),face1) = face1;
        EToE(elem(2),face2) = elem(2);
        EToF(elem(2),face2) = face2;
    end
end

% Algo 1
% $$$ % For each boundary surface element/face
% $$$ for ibc = 1:size(BClines,1)
% $$$
% $$$     % Get boundary vertices (2 vertices in 2D, 3 in 3D with tets)
% $$$     bcnode1 = BClines(ibc,1);
% $$$     bcnode2 = BClines(ibc,2);
% $$$
% $$$     % Find elements sharing same vertices (ie boundary corresponds to element face)
% $$$     elem = [];
% $$$     for ielem = 1:Kel
% $$$         %for jj = 1:Nfaces
% $$$         if any(EToV(ielem,:)==bcnode1) && any(EToV(ielem,:)==bcnode2)
% $$$             elem = [elem ielem];
% $$$         end
% $$$     end
% $$$
% $$$     if length(elem)==2 % open surface, elem(1) and elem(2) on both sides of boundary
% $$$
% $$$         % What if parallel mesh partition? Will we always see elements on
% $$$         % both side on open boundary?
% $$$
% $$$         % Find common face between the two elements (corresponding to open boundary)
% $$$         face1 = find(EToE(elem(1),:)==elem(2)); % face # of element 1
% $$$         face2 = find(EToE(elem(2),:)==elem(1)); % face # of element 2
% $$$
% $$$         % Check that face 'face1' of element 1 connects to face 'face2' of element 2
% $$$         assert(EToF(elem(1),face1)==face2)
% $$$         assert(EToF(elem(2),face2)==face1)
% $$$
% $$$         % Manually remove connectivity
% $$$         EToE(elem(1),face1) = elem(1);
% $$$         EToF(elem(1),face1) = face1;
% $$$         EToE(elem(2),face2) = elem(2);
% $$$         EToF(elem(2),face2) = face2;
% $$$     end
% $$$ end

%%

% Build connectivity maps
[vmapM, mapM, vmapP, mapP, vmapB, mapB] = BuildMaps2D(x,y,Fmask);

vmapM = reshape(vmapM,[Nfp*Nfaces,Kel]);
vmapP = reshape(vmapP,[Nfp*Nfaces,Kel]);
mapP = reshape(mapP,[Nfp*Nfaces,Kel]);
mapM = reshape(mapM,[Nfp*Nfaces,Kel]);

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2D_quad(Norder, r, s);
Drw = (Vanderm*Vr')/(Vanderm*Vanderm'); Dsw = (Vanderm*Vs')/(Vanderm*Vanderm');
