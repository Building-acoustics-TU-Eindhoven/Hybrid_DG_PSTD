function [LIFT] = Lift2Dquad(r,s,Fmask)

% function [LIFT] = Lift2Dquad()
% Purpose  : Compute surface to volume lift term for DG formulation

Globals2D;
Emat = zeros(Np, Nfaces*Nfp);

% face 1 : south
faceR = r(Fmask(:,1));
V1D = Vandermonde1D(Norder, faceR);
massEdge1 = inv(V1D*V1D');
Emat(Fmask(:,1),1:Nfp)         = massEdge1;

% face 2 : east
faceS = s(Fmask(:,2));
V1D = Vandermonde1D(Norder, faceS);
massEdge2 = inv(V1D*V1D');
Emat(Fmask(:,2),Nfp+1:2*Nfp)   = massEdge2;

% face 3 : north
faceR = r(Fmask(:,3));
V1D = Vandermonde1D(Norder, faceR);
massEdge3 = inv(V1D*V1D');
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massEdge3;

% face 4 : west
faceS = s(Fmask(:,4));
V1D = Vandermonde1D(Norder, faceS);
massEdge4 = inv(V1D*V1D');
Emat(Fmask(:,4),3*Nfp+1:4*Nfp) = massEdge4;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = Vanderm*(Vanderm'*Emat);
return
