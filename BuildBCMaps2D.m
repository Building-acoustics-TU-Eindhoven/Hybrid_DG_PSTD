function [mapsBC vmapsBC bmapsBC] = BuildBCMaps2D(BCType,BCcodes,vmapM,mapB)

% function BuildBCMaps2D
% Purpose: Build nodal maps from BCType associated with each boundary codes.
%          The corresponding maps are stored in the cell arrays mapsBC and
%          vmapsBC (this way, the vector BCcodes can contain several
%          boundary codes at once)

Globals2D;

% process BCType
bct    = BCType.';
bnodes = ones(Nfp,1)*bct(:).';
bnodes = bnodes(:);

% find matching boundary codes
mapsBC = cell(size(BCcodes));
vmapsBC = cell(size(BCcodes));
bmapsBC = cell(size(BCcodes));

for ii=1:length(BCcodes)

    % map for variables defined on elements faces (similar to mapB)
    mapsBC{ii} = find(bnodes==BCcodes(ii));

    % map for variables defined on whole elements (similar to vmapB)
    vmapsBC{ii} = vmapM(mapsBC{ii});

    % map for variables defined on boundaries (mapB(bmapsBC{ii})==mapsBC{ii})
    [~,bmapsBC{ii},~] = intersect(mapB,mapsBC{ii});
end

% $$$ function [mapsBC vmapsBC bmapsBC] = BuildBCMaps2D(BCType,BCcodes,vmapM,mapB)
% $$$
% $$$ % function BuildBCMaps2D
% $$$ % Purpose: Build nodal maps from BCType associated with each boundary codes.
% $$$ %          The corresponding maps are stored in the cell arrays mapsBC and
% $$$ %          vmapsBC (this way, the vector BCcodes can contain several
% $$$ %          boundary codes at once)
% $$$
% $$$ Globals2D;
% $$$
% $$$ % process BCType
% $$$ bct    = BCType.';
% $$$ bnodes = ones(Nfp,1)*bct(:).';
% $$$ bnodes = bnodes(:);
% $$$
% $$$ % find matching boundary codes
% $$$ mapsBC = cell(size(BCcodes));
% $$$ vmapsBC = cell(size(BCcodes));
% $$$ bmapsBC = cell(size(BCcodes));
% $$$
% $$$ for ii=1:length(BCcodes)
% $$$     % map for variables defined on elements faces (similar to mapB)
% $$$     mapsBC{ii} = find(bnodes==BCcodes(ii));
% $$$
% $$$     % map for variables defined on whole elements (similar to vmapB)
% $$$     vmapsBC{ii} = vmapM(mapsBC{ii});
% $$$
% $$$     % map for variables defined on boundaries (mapB(bmapsBC{ii})==mapsBC{ii})
% $$$     [~,bmapsBC{ii},~] = intersect(mapB,mapsBC{ii});
% $$$ end
