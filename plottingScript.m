% PS
subplot(211)
surf(x_PS,y_PS,pa_PS,'EdgeColor','none');
%xlabel('$x$ (m)','interpreter','latex')
ylabel('$y$ (m)','interpreter','latex')

plot(recx_PS,recy_PS,'.','color',colors(2,:),'markersize',20)

title('PS')

if ~op_hide_cbar
    cbar_PS = colorbar;
    cbar_PS.Label.String = 'pressure, PS';
end

if strcmp(op_color_max,'max')
    caxis(gather(max(abs(op_plot_transf([pa_PS(:); pa(:)]))).*op_color_range))
elseif strcmp(op_color_max,'const')
    caxis(gather(op_color_range))
else
    caxis(gather(op_plot_transf(op_color_max/(1+sqrt(c0*time))).*op_color_range))
end

colormap(op_colormap)
shading(op_shading)
lighting gouraud
axis tight
view(2)
grid off

% display PMLs
if PMLxb_thickness_PS, line([PMLxb_int_PS PMLxb_int_PS],[PMLyb_ext_PS PMLyt_ext_PS],'color',colors(1,:)); end
if PMLxt_thickness_PS, line([PMLxt_int_PS PMLxt_int_PS],[PMLyb_ext_PS PMLyt_ext_PS],'color',colors(1,:)); end
if PMLyb_thickness_PS, line([PMLxb_ext_PS PMLxt_ext_PS],[PMLyb_int_PS PMLyb_int_PS],'color',colors(1,:)); end
if PMLyt_thickness_PS, line([PMLxb_ext_PS PMLxt_ext_PS],[PMLyt_int_PS PMLyt_int_PS],'color',colors(1,:)); end
% display coupling zone
line([xmin_PS xmax_PS],[ymin_PS ymin_PS],'Color','red','LineWidth',1);

hold off

% DG
subplot(212)
trisurf(PLOTMAT,x(:),y(:),op_plot_transf(pa(:)));
hold all

xlabel('$x$ (m)','interpreter','latex')
ylabel('$y$ (m)','interpreter','latex')

% also plot source and micropones
hold on
%plot(xs,ys,'.','color',colors(1,:),'markersize',20)
plot(recx,recy,'.','color',colors(2,:),'markersize',20)

if ~op_hide_cbar
    cbar_DG = colorbar;
    cbar_DG.Label.String = 'pressure, DG';
end

if strcmp(op_color_max,'max')
    caxis(gather(max(abs(op_plot_transf([pa_PS(:); pa(:)]))).*op_color_range))
elseif strcmp(op_color_max,'const')
    caxis(gather(op_color_range))
else
    caxis(gather(op_plot_transf(op_color_max/(1+sqrt(c0*time))).*op_color_range))
end

colormap(op_colormap)
shading(op_shading)
lighting gouraud
axis tight
view(2)
grid off

title('DG')

% display PMLs
if PMLxb_thickness, line([PMLxb_int PMLxb_int],[PMLyb_ext PMLyt_ext],'color',colors(1,:)); end
if PMLxt_thickness, line([PMLxt_int PMLxt_int],[PMLyb_ext PMLyt_ext],'color',colors(1,:)); end
if PMLyb_thickness, line([PMLxb_ext PMLxt_ext],[PMLyb_int PMLyb_int],'color',colors(1,:)); end
if PMLyt_thickness, line([PMLxb_ext PMLxt_ext],[PMLyt_int PMLyt_int],'color',colors(1,:)); end

% display coupling zone
line([xmin xmax],[ymin_PS ymin_PS]-Nel_CZ*dx_PS,'Color','red','LineWidth',1);
hold off
