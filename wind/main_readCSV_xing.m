clear all;
close all;

wind_file = 'HW_04_RANS.csv';
wind_file = 'canyon_xing/Data_BWL.csv'; % NB & BWL

data = readmatrix(wind_file);

x = 12-data(:,4);
y = data(:,3);
%z = data(:,3);

ux = -data(:,6);
uy = data(:,5);
%uz = data(:,6);

% $$$ duxdx = data(:,7);
% $$$ duxdy = data(:,8);
% $$$ duxdz = data(:,9);
% $$$ duydx = data(:,10);
% $$$ duydy = data(:,11);
% $$$ duydz = data(:,12);
% $$$ duzdx = data(:,13);
% $$$ duzdy = data(:,14);
% $$$ duzdz = data(:,15);

figure
%scatter(x,z,10,sqrt(ux.^2+0*uz.^2))
scatter(x,y,10,ux)
colorbar
axis image
xlabel('$x$ (m)','interpreter','latex')
ylabel('$y$ (m)','interpreter','latex')