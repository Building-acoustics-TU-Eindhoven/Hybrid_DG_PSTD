function [sampleweights,sampleID] = Sample2D(xout, yout)
% function [sampleweights,sampleID] = Sample2D(xout, yout)
% purpose: input = coordinates of output data point
%          output = ID of containing quadrangle and interpolation weights


Globals2D;

% convert xout and yout to column vectors, if need be
if size(xout,1)<size(xout,2)
    xout = xout.';
end
if size(yout,1)<size(yout,2)
    yout = yout.';
end

% coordinates of reference element
v1 = [-1 -1]; v2 = [1 -1]; v3 = [1 1]; v4 = [-1 1];

% find containing element
sampleID = nan(length(xout),1);
pos = nan(length(xout),2); % r,s coordinates in element
for ii=1:Kel
    % split the quadrangle into 2 triangles to facilitate the search
    EToV_local1 = EToV(ii,[1 2 3]);
    EToV_local2 = EToV(ii,[1 3 4]);
    [ID1,B1] = tsearchn([VX.', VY.'], EToV_local1, [xout,yout]);
    [ID2,B2] = tsearchn([VX.', VY.'], EToV_local2, [xout,yout]);

    if any(~isnan([ID1; ID2])) % if one value is not a NaN
        [iout,itri]=find(~isnan([ID1 ID2]));
        sampleID(iout) = ii;

        % barycentric coordinates -> reference element coordinates
        if itri==1
            pos(iout,:) = B1(iout,:)*[v1; v2; v3];
        else
            pos(iout,:) = B2(iout,:)*[v1; v3; v4];
        end
    end
end

rout = pos(:,1); sout = pos(:,2);

% build generalized Vandermonde for the sample point
Vout = Vandermonde2D_quad(Norder, rout, sout);

% build interpolation matrix for the sample point
sampleweights = Vout*invVanderm;

return;
