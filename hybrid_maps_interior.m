% Function for mapping nodes between PSTD and DG
% This code only maps the internal (IN) DG nodes, i.e. only the non-boundary
% nodes of each DG element
% DG QUADS elements

function [map_DG_INxy, map_PS_INxy, x_PS_INxy, y_PS_INxy] = ...
    hybrid_maps_interior(x, y, x_PS, y_PS, dx_PS, xmin, xmax, ymin, ymax)


Globals2D


% Find the inteporlating distances
[r,~] = Nodes2D_quad(Norder); r = r(:);
int_dist_vec = unique((1+r)*dx_PS/2);
N_dist = length(int_dist_vec);

int_dist_vec1 = int_dist_vec(2:end-1);
N_dist1 = length(int_dist_vec1);
tol_1 = min(int_dist_vec1)/10;


% Limits domain
x_min_DG_PStoDG = xmin - tol_1;
x_max_DG_PStoDG = xmax + tol_1;
y_min_DG_PStoDG = ymin - tol_1;
y_max_DG_PStoDG = ymax + tol_1;


% Building interpolating PSTD grids and limits
for iii = 1:N_dist1
    for jjj = 1:N_dist1

        % Interpolation in the XY direction (interpolation in the + direction)
        x_PS_INxy{iii,jjj} = x_PS + int_dist_vec1(jjj);
        y_PS_INxy{iii,jjj} = y_PS + int_dist_vec1(iii);

        % PSTD limits of the WEST, EAST, SOUTH and NORTH copy zones for
        % different interpolating distances
        x_min_PSTD_PStoDG_int_xplus{iii,jjj} = xmin + int_dist_vec1(jjj);
        x_max_PSTD_PStoDG_int_xplus{iii,jjj} = xmax + int_dist_vec1(jjj);

        y_min_PSTD_PStoDG_int_yplus{iii,jjj} = ymin + int_dist_vec1(iii);
        y_max_PSTD_PStoDG_int_yplus{iii,jjj} = ymax + int_dist_vec1(iii);

    end
end


% PSTD indexes
iPS_cozo_PSDG_x_stg_xplus = find(x_min_PSTD_PStoDG_int_xplus{1,1} <= x_PS_INxy{1,1} &...
                                 x_PS_INxy{1,1} <= x_max_PSTD_PStoDG_int_xplus{1,1});
iPS_cozo_PSDG_y_stg_xplus = find(y_min_PSTD_PStoDG_int_yplus{1,1} <= y_PS_INxy{1,1}(iPS_cozo_PSDG_x_stg_xplus) &...
                                 y_PS_INxy{1,1}(iPS_cozo_PSDG_x_stg_xplus) <= y_max_PSTD_PStoDG_int_yplus{1,1});
iPS_cozo_PSDG_stg_xplus = iPS_cozo_PSDG_x_stg_xplus(iPS_cozo_PSDG_y_stg_xplus);

% DG indexes
iDG_cozo_PSDG_x_stg_xplus = find(x_min_DG_PStoDG <= x & x <= x_max_DG_PStoDG);
iDG_cozo_PSDG_y_stg_xplus = find(y_min_DG_PStoDG <= y(iDG_cozo_PSDG_x_stg_xplus) &...
                                 y(iDG_cozo_PSDG_x_stg_xplus) <= y_max_DG_PStoDG);
iDG_cozo_PSDG_stg_xplus = iDG_cozo_PSDG_x_stg_xplus(iDG_cozo_PSDG_y_stg_xplus);


% Finding coincident nodes between PSTD and DG meshes
% Mapping indexes
clear map_DG_PStoDG_INxy map_PS_PStoDG_INxy
for iii = 1:N_dist1
    for jjj = 1:N_dist1
        clear map_DG_PStoDG1_stg_xplus
        for ii = 1:length(iPS_cozo_PSDG_stg_xplus)
            map_DG_PStoDG1_stg_xplus{ii} = find(abs(x(iDG_cozo_PSDG_stg_xplus)-x_PS_INxy{iii,jjj}(iPS_cozo_PSDG_stg_xplus(ii))) < tol_1 &...
                                                abs(y(iDG_cozo_PSDG_stg_xplus)-y_PS_INxy{iii,jjj}(iPS_cozo_PSDG_stg_xplus(ii))) < tol_1);
        end

        [nrows_stg_xplus, ncols_stg_xplus] = cellfun(@size, map_DG_PStoDG1_stg_xplus);
        map_PS_PStoDG1_stg_xplus = cell(1,length(iPS_cozo_PSDG_stg_xplus));
        for ii = 1:length(iPS_cozo_PSDG_stg_xplus)
            map_PS_PStoDG1_stg_xplus(1,ii) = {ones(nrows_stg_xplus(ii),1)*iPS_cozo_PSDG_stg_xplus(ii)};
        end

        map_DG_INxy{iii,jjj} = iDG_cozo_PSDG_stg_xplus(vertcat(map_DG_PStoDG1_stg_xplus{:}));
        map_PS_INxy{iii,jjj} = vertcat(map_PS_PStoDG1_stg_xplus{:});

    end
end

% figure(10)
% % RPM_PlotMesh2D()
% hold on
% for iii = 1:N_dist1
%     for jjj = 1:N_dist1
%         plot(x(map_DG_INxy{iii,jjj}),y(map_DG_INxy{iii,jjj}),'o','Color',[0 0 1],'MarkerEdgeColor',[0 0 1])
%         plot(x_PS_INxy{iii,jjj}(map_PS_INxy{iii,jjj}),y_PS_INxy{iii,jjj}(map_PS_INxy{iii,jjj}),'.')
%         axis equal
%     end
% end

return