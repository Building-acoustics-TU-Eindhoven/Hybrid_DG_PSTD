function BCType = BuildBCType2D(BClines,vmapM,mapB)
% Purpose: Build BCType table, which is a (Kel x Nfaces) matrix containing
% the boundary code associated with each 2D element faces, if any.
% BClines is a matrix associated with each 1D boundary elements
% containing the node indices and the boundary code.

Globals2D;

% $$$ % find elements on the boundaries (just to reduce the number of elements for
% $$$ % the search)
% $$$ [~,BCel] = ind2sub(size(vmapM),mapB);

% initialize BCType matrix (-1 is the default value, corresponding to
% 'not a boundary')
BCType = -1*ones(Kel,Nfaces);

% loop over the 1D boundary elements
for ii=1:size(BClines,1)
    for kk=1:Kel

        % for each BClines, find matching elements
        if sum(ismember(BClines(ii,1:2), EToV(kk,:)))==2

            % find position within element
            id1 = find(EToV(kk,:)==BClines(ii,1));
            id2 = find(EToV(kk,:)==BClines(ii,2));

            % find the face (the ids are always consecutive due to the node reordering)
            % [1 2] -> 1*2 = 2  -> face 1
            % [2 3] -> 2*3 = 6  -> face 2
            % [3 4] -> 3*4 = 12 -> face 3
            % [4 1] -> 4*1 = 4  -> face 4
            switch id1*id2
              case 2
                iface = 1;
              case 6
                iface = 2;
              case 12
                iface = 3;
              case 4
                iface = 4;
            end

            % add boundary code
            BCType(kk,iface) = BClines(ii,3);
            break

        end
    end

end
