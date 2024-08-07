function [mapPatchIntr, mapPatchExtr, mapPatchBndr] = ...
    AcousticsPatchSetup2D(xs, ys, Nlevels)

Globals2D;

ids = [];

% find element containing the source
[~, ids] = Sample2D_quad(xs, ys);


%% grow the patch through 2 neighbors
%Nlevels = 3;
for lev=1:Nlevels
    % find elements connected to element ids
    ids = unique(union(ids, EToE(ids,:)));
end

% not a connectivity matrix, just says if it's in the patch or not
EToPatch = zeros(Kel,1);
EToPatch(ids) = 1;

%% list of surface nodes on
mapPatchIntr = [];
mapPatchExtr = [];
mapPatchBndr = [];
for kM=1:Kel
    for fM=1:Nfaces
        kP = EToE(kM,fM);
        fP = EToF(kM,fM);
        idM = (kM-1)*Nfp*Nfaces+(fM-1)*Nfp + (1:Nfp); % generate unique id for each node on each face of each element

        if(kP~=-1)
        	if(EToPatch(kM)==1 && kM==kP)  % in patch + on a boundary
              mapPatchBndr = [mapPatchBndr; idM];
            elseif(EToPatch(kM)==1 && EToPatch(kP)==0) % current el in patch, adjacent el not in patch
          	  mapPatchIntr = [mapPatchIntr; idM];
        	elseif(EToPatch(kM)==0 && EToPatch(kP)==1) % current el not in patch, adjacent el in patch
          	  mapPatchExtr = [mapPatchExtr; idM];
        	end
        end
    end
end

return
end
