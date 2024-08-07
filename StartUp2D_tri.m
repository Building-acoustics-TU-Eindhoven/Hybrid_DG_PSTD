% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
Nfp = Norder+1; Np = (Norder+1)*(Norder+2)/2; Nfaces=3;

% Compute nodal set
[x,y] = Nodes2D_tri(Norder); [r,s] = xytors_tri(x,y);

% Build reference element matrices
Vanderm = Vandermonde2D_tri(Norder,r,s); invVanderm = inv(Vanderm);
%MassMatrix = invVanderm'*invVanderm;
[Dr,Ds] = Dmatrices2D_tri(Norder, r, s, Vanderm);

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

% find all the nodes that lie on each edge
fmask1   = find( abs(s+1) < NODETOL)';
fmask2   = find( abs(r+s) < NODETOL)';
fmask3   = find( abs(r+1) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3]';
%Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);

% Create surface integral terms
LIFT = Lift2D_tri(r,s,Fmask);

% calculate geometric factors
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);

% calculate geometric factors
[nx, ny, sJ] = Normals2D_tri(x,y,Dr,Ds,Fmask);
Fscale = sJ./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = tiConnect2D(EToV);

% Build connectivity maps
[vmapM, mapM, vmapP, mapP, vmapB, mapB] = BuildMaps2D(x,y,Fmask);

vmapM = reshape(vmapM,[Nfp*Nfaces,Kel]);
vmapP = reshape(vmapP,[Nfp*Nfaces,Kel]);
mapP = reshape(mapP,[Nfp*Nfaces,Kel]);
mapM = reshape(mapM,[Nfp*Nfaces,Kel]);

% Compute weak operators (could be done in preprocessing to save time)
%[Vr, Vs] = GradVandermonde2D_tri(Norder, r, s);
%Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');
