function [QUAD] = PlotField2D_quad(K, Nout, xin, yin)

% function [QUAD] = PlotField2Dquad(K, Nout, xin, yin)
% Purpose: filled contour plot of solution data


 % build equally spaced grid on reference triangle
Npout = (Nout)*(Nout);

EToVlocal = zeros((Nout-1)*(Nout-1),1);
count = 0;
for i = 1 : Nout - 1
    for j = 1 : Nout - 1;
        count = count + 1;
        Gidx = j + (i-1)*Nout;
        EToVlocal(count,[1 2 3 4]) = [Gidx Gidx+Nout Gidx+1+Nout Gidx+1];
    end
end

% build triangulation for all equally spaced nodes on all elements
QUAD = [];
for k=1:K
  QUAD = [QUAD; EToVlocal+(k-1)*Npout];
end

return
