function [nx, ny, sJ] = Normals2D_quad(x,y,Dr,Ds,Fmask)
% function [nx, ny, sJ] = Normals2Dquad()
% Purpose : Compute outward pointing normals at elements faces and surface Jacobians
% By Allan P. Engsig-Karup.

Globals2D;
xr = Dr*x; yr = Dr*y; xs = Ds*x; ys = Ds*y; %J = xr.*ys-xs.*yr;

% interpolate geometric factors to face nodes
fxr = xr(Fmask, :); fxs = xs(Fmask, :); fyr = yr(Fmask, :); fys = ys(Fmask, :);

% build normals
nx = zeros(Nfaces*Nfp, Kel); ny = zeros(Nfaces*Nfp, Kel);
fid1 = (1:Nfp)';         fid2 = (Nfp+1:2*Nfp)';
fid3 = (2*Nfp+1:3*Nfp)'; fid4 = (3*Nfp+1:4*Nfp)';

% face 1 : south
nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);

% face 2 : east
nx(fid2, :) =  fys(fid2, :); ny(fid2, :) = -fxs(fid2, :);

% face 3 : north
nx(fid3, :) = -fyr(fid3, :); ny(fid3, :) =  fxr(fid3, :);

% face 4 : west
nx(fid4, :) = -fys(fid4, :); ny(fid4, :) =  fxs(fid4, :);

% normalise
sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
return;
