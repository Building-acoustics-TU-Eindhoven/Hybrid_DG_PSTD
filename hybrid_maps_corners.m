% Function for mapping nodes between PSTD and DG
% This code only maps the corner (C) DG nodes, i.e. excluding the
% internal and the boundary nodes of each DG element
% DG QUADRILATERAL elements

function [map_DG_C, map_PS_C] = hybrid_maps_corners(x, y, x_PS, y_PS, xmin, xmax, ymin, ymax)

Globals2D

% adjust domain bounds to prevent rounding errors
xmin = xmin - NODETOL;
xmax = xmax + NODETOL;
ymin = ymin - NODETOL;
ymax = ymax + NODETOL;

% find PS points inside prescribed domain
idx_PS = find(x_PS>xmin & x_PS<xmax & y_PS>ymin & y_PS<ymax);

% find DG nodes inside prescribed domain
%idx_DG = find(x>xmin & x<xmax & y>ymin & y<ymax);

% find DG elements inside prescribed domain
centroidx = mean(VX(EToV),2);
centroidy = mean(VY(EToV),2);
k_DG = find(centroidx>xmin & centroidx<xmax & centroidy>ymin & centroidy<ymax);
row = repmat((1:Np).',1,length(k_DG));
col = repmat(k_DG.',Np,1);
idx_DG = sub2ind(size(x),row,col);

% Find coincident nodes between PS and DG
for ii = 1:length(idx_PS)
    map_DG{ii} = find(abs(x(idx_DG)-x_PS(idx_PS(ii))) < NODETOL & ...
                      abs(y(idx_DG)-y_PS(idx_PS(ii))) < NODETOL);
end
[nmatches] = cellfun(@length, map_DG); % number of matches per PS point

% Build corresponding PS map
map_PS = cell(1,length(idx_PS));
for ii = 1:length(idx_PS)
    map_PS(1,ii) = {ones(nmatches(ii),1)*idx_PS(ii)};
end

map_DG_C = idx_DG(vertcat(map_DG{:}));
map_PS_C = vertcat(map_PS{:});

% safety check
if length(map_PS_C)<size(x_PS,2)
    error(['Something is wrong. Check DG and PS grids within coupling ' ...
           'zone and check PS spatial step.'])
end

% $$$ figure(101)
% $$$ % % RPM_PlotMesh2D()
% $$$ % hold on
% $$$ plot(x(map_DG_C),y(map_DG_C),'o','Color',[0 0 0],'MarkerEdgeColor',[0 0 0])
% $$$ hold on
% $$$ plot(x_PS(map_PS_C),y_PS(map_PS_C),'.')
% $$$  axis equal

return