function [x,y] = QuadMap(xc,yc,xi,eta)

% function [x,y] = QuadMap(xc,yc,xi,eta)
% Note: reordered code output in order to simplify the algorithm.
% Implemented by Allan P. Engsig-Karup, apek@imm.dtu.dk.
%
% Example: Plot quadrilateral
%
% xi = linspace(-1,1,5); eta = xi;
% xc = [1 4 6 2]; yc = [0 1 4 3];
% [XI,ETA] = meshgrid(xi,eta);
% [X,Y] = QuadMap(xc,yc,XI,ETA);
% plot(X,Y,'k.')
% axis off

% vertice coordinates of quadrilateral
x1 = xc(1); y1 = yc(1);
x2 = xc(2); y2 = yc(2);
x3 = xc(3); y3 = yc(3);
x4 = xc(4); y4 = yc(4);

% bilinear mapping form reference square to strightsided quadrilateral
h1 = 1-xi; h2 = 1-eta;
h3 = 1+xi; h4 = 1+eta;
x = 0.25*(x1*h1.*h2+x2*h3.*h2+x3*h3.*h4+x4*h1.*h4);
y = 0.25*(y1*h1.*h2+y2*h3.*h2+y3*h3.*h4+y4*h1.*h4);
return
