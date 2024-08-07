% quantity to plot
qqq = sqrt(ux.^2+uy.^2);
if ~op_DG_only
    qqq_PS = sqrt(ux_PS.^2+uy_PS.^2);
end
% $$$ qqq = duxdx;
% $$$ qqq_PS = duxdx_PS;

qqq_colormap = 'parula';

if op_plot_split

    % PS
    subplot(211)
    surf(x_PS,y_PS,qqq_PS,'EdgeColor','none');
    %xlabel('$x$ (m)','interpreter','latex')
    ylabel('$y$ (m)','interpreter','latex')

    cbar = colorbar;


    title('PS')

    colormap(qqq_colormap)
    shading(op_shading)
    lighting gouraud
    axis image
    view(2)
    grid off

    % display coupling zone
    line([xmin_PS xmax_PS],[ymin_PS ymin_PS],'Color',colors(2,:),'LineWidth',2);
    hold off

    % DG
    subplot(212)
    trisurf(PLOTMAT,x(:),y(:),qqq(:));
    xlabel('$x$ (m)','interpreter','latex')
    ylabel('$y$ (m)','interpreter','latex')

    cbar = colorbar;

    title('DG')

    colormap(qqq_colormap)
    shading(op_shading)
    lighting gouraud
    axis image
    view(2)
    grid off

    % display coupling zone
    line([xmin xmax],[ymin_PS ymin_PS]-Nel_CZ*dx_PS,'Color',colors(2,:),'LineWidth',1);
    hold off

else

    % DG
    trisurf(PLOTMAT,x(:),y(:),qqq(:));
    hold all

    if ~op_DG_only
        surf(x_PS(1+Nel_CZ:end,:),y_PS(1+Nel_CZ:end,:),qqq_PS(1+Nel_CZ:end,:),'EdgeColor','none');
    end
    xlabel('$x$ (m)','interpreter','latex')
    ylabel('$y$ (m)','interpreter','latex')

    cbar = colorbar;
    cbar.Label.String = 'wind modulus (m/s)';

    colormap(qqq_colormap)
    shading(op_shading)
    lighting gouraud
    axis image
    view(2)
    grid off

end % op_plot_split
