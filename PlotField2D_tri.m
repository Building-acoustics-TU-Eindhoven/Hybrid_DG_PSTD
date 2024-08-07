function [TRI] = PlotField2D_tri(K, Nout, xin, yin)

% function [TRI] = PlotField2D_tri(Nout, xin, yin)
% Purpose: filled contour plot of solution data


% build equally spaced grid on reference triangle
Npout = (Nout+1)*(Nout+2)/2;
sk = 1;
for n=1:Nout+1
  for m=1:Nout+2-n
    counter(n,m) = sk; sk = sk+1;
  end
end

% build triangulation of equally spaced nodes on reference triangle
tri = [];
for n=1:Nout+1
  for m=1:Nout+1-n,
    v1 = counter(n,m);   v2 = counter(n,m+1);
    v3 = counter(n+1,m); v4 = counter(n+1,m+1);
    if(v4)
      tri = [tri;[[v1 v2 v3];[v2 v4 v3]]];
    else
      tri = [tri;[[v1 v2 v3]]];
    end
  end
end

% build triangulation for all equally spaced nodes on all elements
TRI = [];
for k=1:K
  TRI = [TRI; tri+(k-1)*Npout];
end

return
