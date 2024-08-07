function [dmodedr, dmodeds] = GradSimplex2DP_quad(r,s,id,jd)

% function [dmodedr, dmodeds] = GradSimplex2DP(r,s,id,jd)
% Purpose: Return the derivatives of the modal basis (id,jd) on the 2D simplex at (a,b).

fr = JacobiP(r, 0, 0, id); dfr = GradJacobiP(r, 0, 0, id);
gs = JacobiP(s, 0, 0, jd); dgs = GradJacobiP(s, 0, 0, jd);

% r-derivative
dmodedr = dfr.*gs;

% s-derivative
dmodeds = fr.*dgs;
return;
