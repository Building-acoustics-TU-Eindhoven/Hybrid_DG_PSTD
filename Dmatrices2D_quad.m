function [Dr,Ds] = Dmatrices2D_quad(N,r,s,V)

% function [Dr,Ds] = Dmatrices2Dquad(N,r,s,V)
% Purpose : Initialize the (r,s) differentiation matrices
%	    on the simplex, evaluated at (r,s) at order N

[Vr, Vs] = GradVandermonde2D_quad(N, r, s);
Dr = Vr/V; Ds = Vs/V;
return;
