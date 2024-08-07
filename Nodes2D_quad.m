function [x,y] = Nodes2D_quad(N)
% function [x,y] = Nodes2Dquad(N);
% Purpose  : Compute (x,y) nodes in equilateral quadrilateal for polynomial of order N
% By Allan P. Engsig-Karup.

[x]   = JacobiGL(0,0,N); % Legendre-Gauss-Lobatto nodes
[x,y] = meshgrid(x,x); % tensorproduct grid

return
