function [vmapM, mapM, vmapP, mapP, vmapB, mapB] = BuildMaps2D(x,y,Fmask)

% function [vmapM, mapM, vmapP, mapP, vmapB, mapB] = BuildMaps2D()
% Purpose: Connectivity and boundary tables in the K # of Np elements

Globals2D;

oEToE=EToE;
oEToF=EToF;
oEmap=1:Kel;


% number volume nodes consecutively
nodeids = reshape(1:Kel*Np, Np, Kel);
nodeids = nodeids(:,oEmap);

oKel = size(oEToE,1);

vmapM = zeros(Nfp, Nfaces, oKel);
vmapP = zeros(Nfp, Nfaces, oKel);
mapM = (1:oKel*Nfp*Nfaces).';
mapP = reshape(mapM, Nfp, Nfaces, oKel);

% find index of face nodes with respect to volume node ordering
for k1=1:oKel
    for f1=1:Nfaces
        vmapM(:,f1,k1) = nodeids(Fmask(:,f1), k1);
    end
end

one = ones(1, Nfp);
for k1=1:oKel
    for f1=1:Nfaces

        % find neighbor
        jj = oEToE(k1,f1);
        k2 = find(oEmap==jj); f2 = oEToF(k1,f1);

        % reference length of edge
        v1 = EToV(k1,f1); v2 = EToV(k1, 1+mod(f1,Nfaces));
        refd = sqrt( (VX(v1)-VX(v2))^2 + (VY(v1)-VY(v2))^2 );

        % find volume node numbers of left and right nodes
        vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);
        x1 = x(vidM); y1 = y(vidM); x2 = x(vidP); y2 = y(vidP);
        x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

        % Compute distance matrix
        D = (x1 -x2').^2 + (y1-y2').^2;
        [idM, idP] = find(sqrt(abs(D))<NODETOL*refd);
        vmapP(idM,f1,k1) = vidP(idP); mapP(idM,f1,k1) = idP + (f2-1)*Nfp+(k2-1)*Nfaces*Nfp;
    end
end

% reshape vmapM and vmapP to be vectors and create boundary node list
vmapP = vmapP(:); vmapM = vmapM(:); mapP = mapP(:);
mapB = find(vmapP==vmapM); vmapB = vmapM(mapB);
return
