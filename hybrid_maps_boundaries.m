% Function for mapping nodes between PS and DG
% This code only maps the boundary (B) DG nodes, i.e. excluding the
% internal and the corner nodes of each DG element
% DG QUADRILATERAL elements

function [map_DGx, map_DGy, map_PSx, map_PSy] = hybrid_maps_boundaries(x, y, x_PS, y_PS, dx_PS, xmin, xmax, ymin, ymax)


Globals2D

% Find the interpolating distances
[r,~] = Nodes2D_quad(Norder); r = r(:);
int_dist_vec = unique((1+r)*dx_PS/2);
int_dist_vec1 = int_dist_vec(2:end-1);
N_dist1 = length(int_dist_vec1);

% Limits domain
xmin = xmin - NODETOL;
xmax = xmax + NODETOL;
ymin = ymin - NODETOL;
ymax = ymax + NODETOL;

% DG indexes
%iDG = find(xmin<x & x<xmax & y>ymin & y<ymax);
centroidx = mean(VX(EToV),2);
centroidy = mean(VY(EToV),2);
k_DG = find(centroidx>xmin & centroidx<xmax & centroidy>ymin & centroidy<ymax);
row = repmat((1:Np).',1,length(k_DG));
col = repmat(k_DG.',Np,1);
iDG = sub2ind(size(x),row,col);

% Building interpolating PS grids and limits for the x direction
for iii = 1:N_dist1

    % Interpolation in the x direction (in the + direction)
    x_PSx{iii} = x_PS + int_dist_vec1(iii);
    y_PSx{iii} = y_PS;

    % PS limits of the copy zone for different interpolating distances
    xmin_PSx{iii} = xmin + int_dist_vec1(iii);
    xmax_PSx{iii} = xmax + int_dist_vec1(iii);

    ymin_PSx{iii} = ymin;
    ymax_PSx{iii} = ymax;

end

% PS indexes
iPSx = find(x_PSx{1}>xmin_PSx{1} & x_PSx{1}<xmax_PSx{1} & ...
            y_PSx{1}>ymin_PSx{1} & y_PSx{1}<ymax_PSx{1});

% Finding coincident nodes between PS and DG meshes
% Mapping indexes
clear map_DG_tmp map_PS_tmp
for iii = 1:N_dist1
    clear map_DG_tmp
    for ii = 1:length(iPSx)
        map_DG_tmp{ii} = find(abs(x(iDG)-x_PSx{iii}(iPSx(ii))) < NODETOL &...
                              abs(y(iDG)-y_PSx{iii}(iPSx(ii))) < NODETOL);
    end

    [nrowsx, ncolsx] = cellfun(@size, map_DG_tmp);
    map_PS_tmp = cell(1,length(iPSx));
    for ii = 1:length(iPSx)
        map_PS_tmp(1,ii) = {ones(nrowsx(ii),1)*iPSx(ii)};
    end

    map_DGx{iii} = iDG(vertcat(map_DG_tmp{:}));
    map_PSx{iii} = vertcat(map_PS_tmp{:});

end



% Building interpolating PS grids and limits for the y direction
for iii = 1:N_dist1

    % Interpolation in the y direction (interpolation in the + direction)
    x_PSy{iii} = x_PS;
    y_PSy{iii} = y_PS + int_dist_vec1(iii);

    % PS limits of the copy zone for different interpolating distances
    xmin_PSy{iii} = xmin;
    xmax_PSy{iii} = xmax;

    ymin_PSy{iii} = ymin + int_dist_vec1(iii);
    ymax_PSy{iii} = ymax + int_dist_vec1(iii);

end

% PS indexes
iPSy = find(x_PSy{1}>xmin_PSy{1} & x_PSy{1}<xmax_PSy{1} & ...
            y_PSy{1}>ymin_PSy{1} & y_PSx{1}<ymax_PSy{1});

% Finding coincident nodes between PS and DG meshes
% Mapping indexes
clear map_DG_tmp map_PS_tmp
for iii = 1:N_dist1
    clear map_DG_tmp
    for ii = 1:length(iPSy)
        map_DG_tmp{ii} = find(abs(x(iDG)-x_PSy{iii}(iPSy(ii))) < NODETOL &...
                              abs(y(iDG)-y_PSy{iii}(iPSy(ii))) < NODETOL);
    end

    [nrowsy, ncolsy] = cellfun(@size, map_DG_tmp);
    map_PS_tmp = cell(1,length(iPSy));
    for ii = 1:length(iPSy)
        map_PS_tmp(1,ii) = {ones(nrowsy(ii),1)*iPSy(ii)};
    end

    map_DGy{iii} = iDG(vertcat(map_DG_tmp{:}));
    map_PSy{iii} = vertcat(map_PS_tmp{:});

end


% figure(10)
% % RPM_PlotMesh2D()
% for iii = 1:N_dist1
% hold on
%     plot(x(map_DGx{1,iii}),y(map_DG_Bx{1,iii}),'o','Color',[1 0 0],'MarkerEdgeColor',[1 0 0])
%     plot(x(map_DG_By{1,iii}),y(map_DG_By{1,iii}),'o','Color',[1 0 0],'MarkerEdgeColor',[1 0 0])
%     plot(x_PS_Bx{1,iii}(map_PS_Bx{1,iii}),y_PS_Bx{1,iii}(map_PS_Bx{1,iii}),'.')
%     plot(x_PS_By{1,iii}(map_PS_By{1,iii}),y_PS_By{1,iii}(map_PS_By{1,iii}),'.')
%     axis equal
% end

return