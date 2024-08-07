function [sampleweights,sampleID] = Sample2D_tri(xout, yout)
% function [sampleweights,sampleID] = Sample2D_tri(xout, yout)
% purpose: input = coordinates of output data point
%          output = ID of containing triangle and interpolation weights

Globals2D;

% convert xout and yout to column vectors, if need be
if size(xout,1)<size(xout,2)
    xout = xout.';
end
if size(yout,1)<size(yout,2)
    yout = yout.';
end

% find containing element
[sampleID,B] = tsearchn([VX.', VY.'], EToV, [xout,yout]);

% barycentric coordinates -> reference element coordinates
v1=[-1 -1]; v2=[1 -1]; v3 = [-1 1]; % coordinates of reference element
pos = B*[v1; v2; v3];

rout = pos(:,1); sout = pos(:,2);

% build generalized Vandermonde for the sample point
Vout = Vandermonde2D_tri(Norder, rout, sout);

% build interpolation matrix for the sample point
sampleweights = Vout*invVanderm;

return;
