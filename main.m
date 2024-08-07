% The variables related to DG don't have any suffix to improve
% The variables with the '_PS' suffix are related to PSTD
% readability of the code

clear all;
colors = get(gca,'colororder');
close all;
clc;

set(0,'DefaultFigureWindowStyle','docked')


addpath(genpath(fullfile('.','ext')))
addpath(genpath(fullfile('.','wind')))

Globals2D;

rho0=1.2; % air density
c0=344; % sound speed

totime = 10/c0; % total physical duration (in seconds) (not exact, it's adjusted later based on time step)

%% Simulation settings

op_DG_only = false; % set to true for DG only, otherwise hybrid PSTD-DG model

op_plot = 1; % plot results during execution
op_record = 0; % write pressure signals to disk
op_record_space = 0; % write pressure field over the WHOLE DG domain if
                     % op_snapshot>0 (be careful, in case of large simulations)
op_snapshot = 10; % display (and save/print, if
                  % op_record_space=1/op_print_snapshots=1) pressure field
                  % every n time iterations
op_print_snapshots = 0; % print snapshots as image files if op_snapshot>0
op_plot_mesh = 1; % plot mesh or not (can be costly for large meshes)
op_print_mesh = 0; % print mesh as image file if op_plot_mesh=1

output_name = 'euronoise_test'; % string appended to the output files

% more plot parameters
op_plot_config = 'default'; % select predefined preset
if strcmp(op_plot_config,'seqhot') % sequential colors
    op_plot_split = 0; % split DG and PS plots, or merge them into one plot
    op_colormap = 'hot'; % colormap style
    op_shading = 'interp'; % 'flat' or 'interp'
    op_hide_cbar = 1;
    op_color_max = 0.15; % 'max' or a number
    op_color_range = [0 1]; % [-1 1] for diverging colormap, [0 1] for sequential
    op_plot_transf = @(x) abs(x);
else % default
    op_plot_split = 1; % split DG and PS plots, or merge them into one plot
    op_colormap = 'jet'; % colormap style
    op_shading = 'interp'; % 'flat' or 'interp'
    op_hide_cbar = 0;
    op_color_max = 'max'; % 'max' or a number
    op_color_range = [-1 1]; % [-1 1] for diverging colormap, [0 1] for sequential
    op_plot_transf = @(x) x;
end
if op_DG_only, op_plot_split = 0; end

% function used to print snapshots
op_print_snapshotf = @(n,output_name) export_fig( ...
    fullfile('.','outputs','snapshots',[output_name '_' num2str(n) '.png']),'-m2.5','-silent');

% DG parameters
Norder = 4; % polynomial interpolation order for DG
Nout = Norder; % interpolation order for plots

% DG mesh file
msh_fname = 'euronoise_wing';
%gmsh_opts = '-setnumber lc .5';

% wind input data
%wind_fname = 'HW_04_RANS.csv';
wind_fname = '';

% run GMSH from Matlab (uncomment if you don't want this, and if the msh already exists)
geo_file = fullfile(pwd,'meshes',[msh_fname '.geo']);
if exist(geo_file,'file')
    disp('Running GMSH')
    gmsh_opts = '';
    system(['gmsh -2 ' gmsh_opts ' ' geo_file]);
else
    error('geo file does not exist!')
end

if ~op_DG_only
    % PS parameters
    dx_PS = 1.5;
    dy_PS = dx_PS;
    fmax_PS = c0/(dx_PS*2);
    lambda_min = c0/fmax_PS;

    % PS geometry
    Ny_PS_above_DG = 101; % number of PS points on top of DG in the y direction
    % Coupling zone (CZ)
    Nel_CZ_PS2DG = 2;
    Nel_CZ_DG2PS = 30;
    % + check mesh structured inside CZ
end

% source position (xs,ys) in meters
xs=45;
ys=12;

% initial conditions (source)
S_0=1; % amplitude
FWHM = 5*dx_PS; % Full Width at Half Maximum

% definition of initial spatial distribution
ini_distrib = @(x,y) S_0*exp(-4*log(2)*( ((x-xs).^2 + (y-ys).^2)/FWHM^2));

% recording locations (in meters)
% DG
recx=[10];
recy=[5]; % line of sight

if ~op_DG_only
    % PS (make sure that the locations correspond to PS grid points)
    recx_PS=[];
    recy_PS=[];
end

%% Read mesh file created with gmsh (accepted msh format: ascii version 2)
msh_file = fullfile(pwd,'meshes',[msh_fname '.msh']);
if ~exist(msh_file,'file')
    error('mesh file does not exist!')
end
disp('Reading mesh file')
[VX,VY,Nv,Kel,EToV,BClines,Ekind] = RPM_meshread(msh_file);
%%
if op_plot_mesh % plot mesh and boundaries
    figure(1)
    for ii=1:Kel % for each element
        tmp = EToV(ii,:); % find node numbers
        tmp = [tmp tmp(1)]; % define closed path for plotting
        plot(VX(tmp),VY(tmp),'k')
        hold all
    end
    xlabel('x [m]')
    ylabel('y [m]')
    axis equal
end

%% initialize DG solver
disp('Initializing DG solver')
if strcmp(Ekind,'triangles'), StartUp2D_tri; end
if strcmp(Ekind,'quads'), StartUp2D_quad; end

% Build modal filtering matrix
Nc = 3/4*Norder; % cutoff mode (2/3*Norder -> impact only upper third of modes)
sp = 64; % filter order
if strcmp(Ekind,'triangles'), [Filtmod] = Filter2D_tri(Nc,sp); end
if strcmp(Ekind,'quads'), [Filtmod] = Filter2D_quad(Nc,sp); end

% process boundary points
disp('Processing boundaries')
BCType = BuildBCType2D(BClines,vmapM,mapB); % build face type matrix from mesh file
BCcodes = unique(BClines(:,3)); % get mesh boundary codes
[mapsBC vmapsBC bmapsBC] = BuildBCMaps2D(BCType,BCcodes,vmapM,mapB); % build boundary maps

% physical boundary codes
BCparams = cell(size(BCcodes));
for ii=1:length(BCcodes)
    % define 'Rinf', 'poles' and 'residues' associated with each boundary code

    % The pole residue model of the reflection coefficient is of the form:
    % R(omega) = Rinf + sum( residues_real/(s-poles_real) ) ...
    %    + sum( residues_cplx/(s-poles_cplx) + conj(residues_cplx)/(s-conj(poles_cplx)));
    % Note: here, by convention, the real part of the poles is
    % negative. For instance:
    % Rinf = 0.5;
    % residues = [1000; 100; 1000;];
    % poles = [-1000; -200+1i*400*pi; -300+1i*1000*pi;];

    switch BCcodes(ii)
      case 1 % rigid / perfectly reflecting BC
        Rinf = 1;
        residues = [];
        poles = [];
        BCparams{ii} = BCparameters(BCcodes(ii),Rinf,poles,residues);

      case 2 % impedance
        % load impedance data
        TDIBC_repo = 'varpor_25_50_1_2250';
        TDIBC_path = fullfile('.','impedance',TDIBC_repo);
        TDIBC.R0 = dlmread(fullfile(TDIBC_path,'D0.inp'));
        TDIBC.residues = dlmread(fullfile(TDIBC_path,'residues.inp'));
        TDIBC.poles = dlmread(fullfile(TDIBC_path,'poles.inp'));

        Rinf = TDIBC.R0;
        residues = TDIBC.residues(:,1);
        poles = TDIBC.poles(:,1);
        BCparams{ii} = BCparameters(BCcodes(ii),Rinf,poles,residues);

      otherwise % nonreflecting BC
        Rinf = 0;
        residues = [];
        poles = [];
        BCparams{ii} = BCparameters(BCcodes(ii),Rinf,poles,residues);
    end
end

% process pole-residue models
[BCcodesInst,BCcodesReal,BCcodesCplx]=BCparameters.sort(BCcodes,BCparams);

% initialize auxiliaxy variables
phi = cell(size(BCparams));
psir = cell(size(BCparams));
psii = cell(size(BCparams));
resphi = cell(size(BCparams));
respsir = cell(size(BCparams));
respsii = cell(size(BCparams));
for ibc=BCcodesReal
    phi{ibc} = zeros(length(mapsBC{ibc}),length(BCparams{ibc}.lambda));
    resphi{ibc} = zeros(size(phi{ibc}));
end
for ibc=BCcodesCplx
    psir{ibc} = zeros(length(mapsBC{ibc}),length(BCparams{ibc}.alpha));
    psii{ibc} = zeros(length(mapsBC{ibc}),length(BCparams{ibc}.alpha));
    respsir{ibc} = zeros(size(psir{ibc}));
    respsii{ibc} = zeros(size(psii{ibc}));
end

if op_plot_mesh % plot boundaries associated with each boundary code
    figure(1)
    hold all
    for ii=1:length(BCcodes)
        plot(x(vmapsBC{ii}),y(vmapsBC{ii}),'.')
        hold all
    end
    % % also plot interpolation nodes
    % figure(1)
    % hold all
    % plot(x,y,'.','MarkerSize',3)
end

%% create additional maps based on mesh layout, for plots
if strcmp(Ekind,'triangles'), [PLOTMAT] = PlotField2D_tri(Kel, Nout, x, y); end
if strcmp(Ekind,'quads'), [PLOTMAT] = PlotField2D_quad(Kel, Nout+1, x, y); end

xmin = min(x(:));
xmax = max(x(:));
ymin = min(y(:));
ymax = max(y(:));

if ~op_DG_only
    xmin_PS = xmin;
    xmax_PS = xmax;
    ymin_PS = ymax;
    ymax_PS = ymin_PS + Ny_PS_above_DG*dy_PS;

    % check that PS and DG grids match
    assert( abs(xmin-xmin_PS)<NODETOL ...
            & abs(xmax-xmax_PS)<NODETOL ...
            & abs(ymax-ymin_PS)<NODETOL, ...
            'PS and DG grids do not match!' )
end

%% Build PS grid
if ~op_DG_only
    if (Nel_CZ_PS2DG+Nel_CZ_DG2PS)*dy_PS>ymin_PS-ymin % coupling zone
                                                      % extends below DG domain
        error('coupling zone extends below DG domain; adjust size of CZ')
        %Nel_CZ_DG2PS = ceil((ymin_PS-ymin)/dy_PS)-Nel_CZ_PS2DG;
    end

    Nel_CZ = Nel_CZ_PS2DG + Nel_CZ_DG2PS; % Number of vertical points in the coupling zone

    xvec_PS = xmin_PS:dx_PS:xmax_PS;
    yvec_PS = ymin_PS-Nel_CZ*dx_PS:dy_PS:ymax_PS;

    [x_PS,y_PS] = meshgrid(xvec_PS,yvec_PS);

    Nx_PS = length(xvec_PS);
    Ny_PS = length(yvec_PS);

    if mod(Nx_PS,2)~=0
        warning(['x dimension of PS grid not a multiple of 2, computations might be slow'])
    end

    if mod(Ny_PS,2)~=0
        warning(['y dimension of PS grid not a multiple of 2, computations might be slow'])
    end

    % wavenumber vectors
    kxvec_PS = wavenumber_vector(Nx_PS, dx_PS);
    kyvec_PS = wavenumber_vector(Ny_PS, dy_PS);

    % corresponding matrices
    kx_PS = repmat(kxvec_PS,Ny_PS,1);
    ky_PS = repmat(kyvec_PS,Nx_PS,1).';
end

%% sound speed matrix
disp('Processing atmospheric parameters')

% sound speed
c_eff = ones(Np,Kel)*c0; % + b_wind*log( (y/y_0) + 1);
c_eff_face = c_eff(vmapM);

if ~op_DG_only
    c_eff_PS = ones(size(x_PS))*c0; % + b_wind*log( (y_PS/y_0) + 1);
end

%% wind components for DG and PS

% initialization to 0
ux = 0*c0*ones(Np,Kel);
uy = 0*c0*ones(Np,Kel);
if ~op_DG_only
    ux_PS = 0*c0*ones(size(x_PS));
    uy_PS = 0*c0*ones(size(x_PS));
end
% wind derivatives
duxdx = 0*c0*ones(Np,Kel);
duxdy = 0*c0*ones(Np,Kel);
duydx = 0*c0*ones(Np,Kel);
duydy = 0*c0*ones(Np,Kel);
if ~op_DG_only
    duxdx_PS = 0*c0*ones(size(x_PS));
    duxdy_PS = 0*c0*ones(size(x_PS));
    duydx_PS = 0*c0*ones(size(x_PS));
    duydy_PS = 0*c0*ones(size(x_PS));
end

% load wind data from file
wind_fpath = fullfile('wind', wind_fname);

if exist(wind_fpath,'file')==2
    wind.data = readmatrix(wind_fpath);

    wind.x = wind.data(:,2)-0.5;
    wind.y = wind.data(:,3);
    wind.ux = wind.data(:,4);
    wind.uy = wind.data(:,5);
    % wind.duxdx = wind.data(:,7);
    % wind.duxdz = wind.data(:,9);
    % wind.duzdx = wind.data(:,13);
    % wind.duzdz = wind.data(:,15);

    % build interpolant
    InterpMethod = 'linear';
    ExtrapMethod = 'nearest';
    wind.ux_int = scatteredInterpolant(wind.x,wind.y,wind.ux,InterpMethod,ExtrapMethod);
    wind.uy_int = scatteredInterpolant(wind.x,wind.y,wind.uy,InterpMethod,ExtrapMethod);
    % wind.duxdx_int = scatteredInterpolant(wind.x,wind.z,wind.duxdx,InterpMethod,ExtrapMethod);
    % wind.duxdz_int = scatteredInterpolant(wind.x,wind.z,wind.duxdz,InterpMethod,ExtrapMethod);
    % wind.duzdx_int = scatteredInterpolant(wind.x,wind.z,wind.duzdx,InterpMethod,ExtrapMethod);
    % wind.duzdz_int = scatteredInterpolant(wind.x,wind.z,wind.duzdz,InterpMethod,ExtrapMethod);

    % plot CFD wind
    figure
    scatter(wind.x,wind.y,10,wind.ux)
    axis image
    xlabel('$x$ (m)','interpreter','latex')
    ylabel('$y$ (m)','interpreter','latex')
    colormap parula
    cbar = colorbar;
    cbar.Label.String = 'raw horizontal wind (m/s)';
    op_print_snapshotf(0,[output_name '_wind']);

    for ii=1:length(x(:))
        ux(ii) = 0*wind.ux_int(x(ii),y(ii));
        uy(ii) = 0*wind.uy_int(x(ii),y(ii));
        % duxdx(ii) = wind.duxdx_int(x(ii),y(ii));
        % duxdy(ii) = wind.duxdz_int(x(ii),y(ii));
        % duydx(ii) = wind.duzdx_int(x(ii),y(ii));
        % duydy(ii) = wind.duzdz_int(x(ii),y(ii));
    end

    if ~op_DG_only
        for ii=1:length(x_PS(:))
            ux_PS(ii) = 0*wind.ux_int(x_PS(ii),y_PS(ii));
            uy_PS(ii) = 0*wind.uy_int(x_PS(ii),y_PS(ii));
            % duxdx_PS(ii) = wind.duxdx_int(x_PS(ii),y_PS(ii));
            % duxdy_PS(ii) = wind.duxdz_int(x_PS(ii),y_PS(ii));
            % duydx_PS(ii) = wind.duzdx_int(x_PS(ii),y_PS(ii));
            % duydy_PS(ii) = wind.duzdz_int(x_PS(ii),y_PS(ii));
        end
    end

    % generate wind from Monin-Obukhov similarity theory (MOS)
    % ustar = 0.5; % frictional velocity (Meppen, [0.03; 0.3])
    % Fb = 0D-3; % buoyancy flux (Meppen, [-0.01; 0.04]
    % z0 = 0.005; % roughness length (Meppen, z0 = 0.005)
    % T0 = 273.15 + 15; % reference surface temperature
    % zref = z0; % zref : reference height
    % Tref = 273.15 + 15; % temperature at zref
    % windphi = 0; % wind horizontal direction w.r.t. the +x direction

    % windx = @(z) mos(z, ustar, Fb, z0, T0, zref, Tref, windphi);
    % dwindxdz = @(z) dmosdz(z, ustar, Fb, z0, T0, zref, Tref, windphi);

    % ux = windx(y);
    % uy = 0*c0*ones(Np,Kel);
    % if ~op_DG_only
    %     ux_PS = windx(y_PS);
    %     uy_PS = 0*c0*ones(size(x_PS));
    % end

    % % wind derivatives
    % duxdx = 0*c0*ones(Np,Kel);
    % duxdy = dwindxdz(y);
    % duydx = 0*c0*ones(Np,Kel);
    % duydy = 0*c0*ones(Np,Kel);
    % if ~op_DG_only
    %     duxdx_PS = 0*c0*ones(size(x_PS));
    %     duxdy_PS = dwindxdz(y_PS);
    %     duydx_PS = 0*c0*ones(size(x_PS));
    %     duydy_PS = 0*c0*ones(size(x_PS));
    % end

    % plot modulus of interpolated wind used in the simulation
    figure
    plottingScript_wind;
    op_print_snapshotf(1,[output_name '_wind']);
    %export_fig(fullfile('.','outputs','snapshots',[output_name '_wind.png']),'-m2.5','-silent');

end


% find maximum wind modulus, to adjust time step
wind_max = max(sqrt(ux(:).^2 + uy(:).^2));
if ~op_DG_only
    wind_max_PS = max(sqrt(ux_PS(:).^2 + uy_PS(:).^2));
    wind_max = max(wind_max,wind_max_PS);
end

% set PML coefficient beta, to improve stability with wind
if ~op_DG_only
    um = mean([ux(:); ux_PS(:)])/max([c_eff(:); c_eff_PS(:)]); % mean x Mach number
else
    um = mean(ux(:))/max(c_eff(:)); % mean x Mach number
end
PMLx_beta = um/(1-um.^2)/c0;


%% Time parameters

precision_DG = 3; % Toulorge thesis. Error magnitude: (1) Emag=1dB; (2) Emag=0.1dB; (3) Emag=0.01dB; (4) Emag=0.001dB; (5) Emag=0.0001dB;
CFL_TOU = CFL_DG_values(precision_DG,Norder);
if strcmp(Ekind,'triangles'), dxscale = min(dtscale2D_tri); end
if strcmp(Ekind,'quads'), dxscale = min(dtscale2D_quad)/3; end
dt_TOU = dxscale*CFL_TOU/(c0+wind_max);
if ~op_DG_only
    CFL_PS = .5;
    dt_PS = dx_PS*CFL_PS/(c0+wind_max);
    dt_factor = ceil(dt_PS/dt_TOU);
    dt = dt_PS/dt_factor;
    nmax = ceil(totime/dt); % total number of time steps
    time = 0; % running time
    nmax_PS = ceil(totime/dt_PS);
    time_PS = 0; % running time
    n = 0; % time iteration counter
    n_PS = 0;
else
    dt = dt_TOU;
    nmax = ceil(totime/dt); % total number of time steps
    time = 0; % running time
    dt_PS = dt; % to simplify the code for DG only
    nmax_PS = nmax;
    dt_factor = 1;
    n_PS = 0;
end

%% Initialize acoustics variables

% DG primitive variables
pa = ini_distrib(x,y);
vx = zeros(Np,Kel);
vy = zeros(Np,Kel);
% DG time integration variables
respa = zeros(Np,Kel);
resvx = zeros(Np,Kel);
resvy = zeros(Np,Kel);
% DG boundary conditions
qcinc = zeros(size(mapB));
qcref = zeros(size(mapB));

if ~op_DG_only
    % PS primitive variables
    pa_PS = ini_distrib(x_PS,y_PS);
    vx_PS = zeros(size(x_PS));
    vy_PS = zeros(size(x_PS));
    % PS time integration variables
    pa0_PS = zeros(size(x_PS));
    vx0_PS = zeros(size(x_PS));
    vy0_PS = zeros(size(x_PS));
end

% find max initial pressure, to check stability
if ~op_DG_only
    pa_max_ini = max(abs([pa(:); pa_PS(:)]));
else
    pa_max_ini = max(abs(pa(:)));
end

%% PML settings
disp('Processing PMLs')

% definition of profiles for DG and PS (here, polynomial)
profile_sigma = @(pos,amp,pwr,thick) amp*(pos/thick).^pwr;
profile_stretch = @(pos,amp,pwr,thick) (1 - (pos/thick).^pwr ); % for buffer zone test


if ~op_DG_only
    % DG layer thickness (in meters)
    PMLxb_thickness = 20*dx_PS; %LEFT PML (x-axis)
    PMLxt_thickness = 20*dx_PS; %RIGHT PML (x-axis)
    PMLyb_thickness = 0*dy_PS; %BOTTOM PML (y-axis)
    PMLyt_thickness = 0; %TOP PML (y-axis)

    % PS layer thickness (in meters)
    PMLxb_thickness_PS = PMLxb_thickness; %LEFT PML (x-axis)
    PMLxt_thickness_PS = PMLxt_thickness; %RIGHT PML (x-axis)
    PMLyb_thickness_PS = 0; %BOTTOM PML (y-axis)
    PMLyt_thickness_PS = 20*dy_PS; %TOP PML (y-axis)

    % DG maximum absorption
    PMLxb_sigma_max = 1/dt_PS;
    PMLxt_sigma_max = 1/dt_PS;
    PMLyb_sigma_max = 1/dt_PS;
    PMLyt_sigma_max = 1/dt_PS;

    % PS maximum absorption
    PMLxb_sigma_max_PS = PMLxb_sigma_max;
    PMLxt_sigma_max_PS = PMLxt_sigma_max;
    PMLyb_sigma_max_PS = 1/dt_PS;
    PMLyt_sigma_max_PS = 1/dt_PS;

else
    % DG layer thickness (in meters)
    PMLxb_thickness=21*dx_PS; %LEFT PML (x-axis)
    PMLxt_thickness=21*dx_PS; %RIGHT PML (x-axis)
    PMLyb_thickness=21*dx_PS; %BOTTOM PML (y-axis)
    PMLyt_thickness=21*dx_PS; %TOP PML (y-axis)

    % DG maximum absorption
    PMLxb_sigma_max = 0.5/dt;
    PMLxt_sigma_max = 0.5/dt;
    PMLyb_sigma_max = 0.5/dt;
    PMLyt_sigma_max = 0.5/dt;
end

% DG polynomial order
PMLxb_pwr = 2;
PMLxt_pwr = 2;
PMLyb_pwr = 2;
PMLyt_pwr = 2;

if ~op_DG_only
    % PS polynomial order
    PMLxb_pwr_PS = PMLxb_pwr;
    PMLxt_pwr_PS = PMLxt_pwr;
    PMLyb_pwr_PS = 2;
    PMLyt_pwr_PS = 2;
end

if ~op_DG_only
    % PS stretching coefficient (note: there is no stretching for DG!)
    PMLxb_stretch_max_PS = 1; % don't change for hybrid method
    PMLxt_stretch_max_PS = 1; % don't change for hybrid method
    PMLyb_stretch_max_PS = 1;
    PMLyt_stretch_max_PS = 2;
end

BuildPMLs2D;
if ~op_DG_only, BuildPMLs2D_PS; end

if op_plot_mesh % plot PMLs
    figure(1)
    hold on
    if PMLxb_thickness, line([PMLxb_int PMLxb_int],[PMLyb_ext PMLyt_ext],'color',colors(2,:),'linewidth',1.5); end
    if PMLxt_thickness, line([PMLxt_int PMLxt_int],[PMLyb_ext PMLyt_ext],'color',colors(2,:),'linewidth',1.5); end
    if PMLyb_thickness, line([PMLxb_ext PMLxt_ext],[PMLyb_int PMLyb_int],'color',colors(2,:),'linewidth',1.5); end
    if PMLyt_thickness, line([PMLxb_ext PMLxt_ext],[PMLyt_int PMLyt_int],'color',colors(2,:),'linewidth',1.5); end
end


%% Check numerical stability of PMLs and TDIBCs

% For TDIBCs, check that all poles are within stability region
for ii=1:length(BCparams)
    xq = real(BCparams{ii}.poles);
    yq = imag(BCparams{ii}.poles);

    % compute amplification factors for DG (no TDIBC in PS)
    g_DG = 1;
    for ii=1:length(rkf84g)
        g_DG = g_DG + rkf84g(ii)*((xq+1i*yq)*dt).^ii;
    end

    if any(abs(g_DG)>1)
        error('TDIBCs will be unstable with current time step.')
    end
end

% PMLs
xq = -max([PMLx_sigma(:); PMLy_sigma(:)]);
if ~op_DG_only
    xq_PS = -max([PMLx_sigma_PS(:); PMLy_sigma_PS(:)]);
end

% compute PML amplification factors
g_DG = 1;
for ii=1:length(rkf84g)
    g_DG = g_DG + rkf84g(ii)*(xq*dt).^ii;
end

if ~op_DG_only
    g_PS = 1;
    for ii=1:length(RKgamma_PS)
        g_PS = g_PS + RKgamma_PS(ii)*(xq_PS*dt_PS).^ii;
    end
    if any(abs(g_DG)>1) || any(abs(g_PS)>1)
        error('PMLs will be unstable with current time step.')
    end
else
    if abs(g_DG)>1
        error('PMLs will be unstable with current time step.')
    end
end


%% Build coupling zone (CZ)

if ~op_DG_only
    % Find interpolating distances
    int_dist_vec = unique((1+r)*dx_PS/2);
    N_dist = length(int_dist_vec);
    int_dist_vec1 = int_dist_vec(2:end-1);
    N_dist1 = length(int_dist_vec1);

    % CZ PS to DG
    xmin_PS2DG = xmin_PS;
    xmax_PS2DG = xmax_PS;
    ymin_PS2DG = ymin_PS - Nel_CZ_PS2DG*dy_PS;
    ymax_PS2DG = ymin_PS;

    % build corresponding maps
    [map_DG_PS2DG_C,map_PS_PS2DG_C] = ...
        hybrid_maps_corners(x,y,x_PS,y_PS,xmin_PS2DG,xmax_PS2DG,ymin_PS2DG,ymax_PS2DG);
    [map_DG_PS2DG_Bx, map_DG_PS2DG_By, map_PS_PS2DG_Bx, map_PS_PS2DG_By] = ...
        hybrid_maps_boundaries(x,y,x_PS,y_PS,dx_PS,xmin_PS2DG, xmax_PS2DG, ymin_PS2DG, ymax_PS2DG);
    [map_DG_PS2DG_Ixy, map_PS_PS2DG_Ixy] =...
        hybrid_maps_interior(x, y, x_PS, y_PS, dx_PS, xmin_PS2DG, xmax_PS2DG, ymin_PS2DG, ymax_PS2DG);

    % CZ DG to PS
    xmin_DG2PS = xmin_PS;
    xmax_DG2PS = xmax_PS;
    ymin_DG2PS = ymin_PS - Nel_CZ*dy_PS;
    ymax_DG2PS = ymin_PS - (Nel_CZ_PS2DG+1)*dy_PS;

    % build corresponding maps
    [map_DG_DG2PS,map_PS_DG2PS] = hybrid_maps_corners(x,y,x_PS,y_PS,xmin_DG2PS,xmax_DG2PS,ymin_DG2PS,ymax_DG2PS);
end


%% Build CZ spatial window

if ~op_DG_only
    alpha_win = (Nel_CZ_DG2PS+3)/14;
    beta_win = 3;
    N_win = Nel_CZ_DG2PS;
    l_win = 0:N_win-1;
    win_y = exp(-alpha_win*log(10).*((l_win-N_win)/N_win).^(2*beta_win));
    win_yvec = ones(size(yvec_PS));
    win_yvec(1:N_win) = win_y;
    win_PS = repmat(win_yvec.',1,length(xvec_PS));
end

%% Build spectral filter

if ~op_DG_only
    % filter settings (see JCP paper)
    ppw_filter_PS = 2.3;
    alpha_filter_PS = 4;
    beta_filter_PS = 4;
    kmax = pi/dx_PS; % Nyquist–Shannon wavenumber
    kc = 2*pi/(ppw_filter_PS*dx_PS); % cut-off wavenumber

    filter_xy = spectral_filter_2D(kx_PS, ky_PS, ppw_filter_PS, alpha_filter_PS, beta_filter_PS, kmax, kc);
    filter_x = spectral_filter_1D(kx_PS, ppw_filter_PS, alpha_filter_PS, beta_filter_PS, kmax, kc);
    filter_y = spectral_filter_1D(ky_PS, ppw_filter_PS, alpha_filter_PS, beta_filter_PS, kmax, kc);
    kx_PS = filter_x.*kx_PS;
    ky_PS = filter_y.*ky_PS;
end

%% precompute some quantities to save computational time
un = nx.*ux(vmapM) + ny.*uy(vmapM); % normal wind on shared intersections
cnxnx=c_eff_face.*nx.*nx;
cnyny=c_eff_face.*ny.*ny;
cnxny=c_eff_face.*nx.*ny;
nxrho = nx/rho0;
nyrho = ny/rho0;
rccnx=rho0*c_eff_face.^2.*nx;
rccny=rho0*c_eff_face.^2.*ny;
rcnxun = rho0*c_eff_face.*nx.*un;
rcnyun = rho0*c_eff_face.*ny.*un;
rcc = rho0*c_eff.^2;
if ~op_DG_only
    rcc_PS = rho0*c_eff_PS.^2;
end

%% Recording

if op_record
    op_record_DG = 1;
    op_record_PS = 1;
else
    op_record_DG = 0;
    op_record_PS = 0;
end
if op_DG_only, op_record_PS = 0; end

% DG

% remove recording locations outside the computational domain
if exist('recx', 'var')==0, recx = []; end
if exist('recy', 'var')==0, recy = []; end
recx_tmp=recx; recy_tmp=recy;
for ii=1:length(recx)
    if (recx(ii)>xmax || recx(ii)<xmin) || ...
            (recy(ii)>ymax || recy(ii)<ymin)
        recx_tmp(ii) = NaN; recy_tmp(ii) = NaN;
    end
end
recx=recx_tmp(~isnan(recx_tmp)); recy=recy_tmp(~isnan(recy_tmp));
if ~length(recx) op_record_DG=0; end

% find in which elements the recording coordinates are found
if op_record_DG
    if strcmp(Ekind,'triangles'), [rec_interp, rec_id] = Sample2D_tri(recx, recy); end
    if strcmp(Ekind,'quads'), [rec_interp, rec_id] = Sample2D_quad(recx, recy); end
    % rec_interp is the interpolating matrix you can multiply times the solution of the
    % acoustics variables of the rec_id element and record in time (t):
    pa_rec=zeros(length(rec_id),nmax_PS);
end

% PS
if ~op_DG_only
    % find recording coordinates
    recx_tmp=recx_PS; recy_tmp=recy_PS;
    rec_idx_PS = []; rec_idy_PS = [];
    for ii=1:length(recx_PS)
        idx = find(xvec_PS<recx_PS(ii)+NODETOL & xvec_PS>recx_PS(ii)-NODETOL);
        idy = find(yvec_PS<recy_PS(ii)+NODETOL & yvec_PS>recy_PS(ii)-NODETOL);
        if length(idx)==1 && length(idy)==1
            rec_idx_PS = [rec_idx_PS idx];
            rec_idy_PS = [rec_idy_PS idy];
            recx_tmp(ii) = xvec_PS(idx);
            recy_tmp(ii) = yvec_PS(idy);
        else % remove locations that do not correspond to a PS grid point
            recx_tmp(ii) = NaN; recy_tmp(ii) = NaN;
        end
    end
    recx_PS=recx_tmp(~isnan(recx_tmp)); recy_PS=recy_tmp(~isnan(recx_tmp));
    if ~length(recx_PS) op_record_PS=0; end

    if op_record_PS
        pa_rec_PS=zeros(length(rec_idx_PS),nmax_PS);
    end

    if ~op_record_DG && ~op_record_PS, op_record=0; end

else
    if ~op_record_DG, op_record=0; end
end

%%
if op_plot_mesh % also plot source and micropones
    % figure(1)
    % hold on
    % plot(xs,ys,'.','color',colors(1,:),'markersize',20)
    if op_print_mesh
        export_fig(fullfile('.','outputs','snapshots',[output_name '_mesh.png']),'-m2.5','-silent');
    end
end

%% GPU computations
% disable GPU computations if the machine can't run them
try
    gpuArray;
    disp('GPU computations enabled')
catch
    gpuArray = @(x) x;
    gather = @(x) x;
    warning('GPU computations disabled')
end

% send all workspace variables of class 'double' to GPU
wvars = whos;
for ii = 1:length(wvars)
    if strcmp(wvars(ii).class,'double')
        name = wvars(ii).name;
        assignin('base', name, gpuArray(evalin('base', name)));
    end
end

%% Initialize plot
if op_plot
    figure(666)
    plottingScript;

    drawnow
    if op_print_snapshots
        op_print_snapshotf(n_PS,output_name);
    end
end % op_plot

%% Run simulation
disp('Starting simulation')
tic
for n_PS=1:nmax_PS % PS time loop

    % display current time iteration
    disp(['time iteration ' num2str(n_PS) '/' num2str(nmax_PS)])

    if ~op_DG_only

        % window PS variables inside coupling zone
        pa_PS = win_PS .* pa_PS;
        vx_PS = win_PS .* vx_PS;
        vy_PS = win_PS .* vy_PS;

        % store acoustic variables for PS time integration
        time0_PS = time_PS;
        pa0_PS = pa_PS;
        vx0_PS = vx_PS;
        vy0_PS = vy_PS;
        PMLx_pa0_PS = PMLx_pa_PS;
        PMLx_vx0_PS = PMLx_vx_PS;
        PMLx_vy0_PS = PMLx_vy_PS;
        PMLy_pa0_PS = PMLy_pa_PS;
        PMLy_vx0_PS = PMLy_vx_PS;
        PMLy_vy0_PS = PMLy_vy_PS;

        for rkstep=1:RKnstage_PS % PS time integration loop

            % compute spatial x-derivatives
            dpadx = PMLx_stretch_PS .* ifft(1i*kx_PS.*fft(pa_PS,[],2),[],2,'symmetric');
            dvxdx = PMLx_stretch_PS .* ifft(1i*kx_PS.*fft(vx_PS,[],2),[],2,'symmetric');
            dvydx = PMLx_stretch_PS .* ifft(1i*kx_PS.*fft(vy_PS,[],2),[],2,'symmetric');

            % compute spatial y-derivatives
            dpady = PMLy_stretch_PS .* ifft(1i*ky_PS.*fft(pa_PS,[],1),[],1,'symmetric');
            dvxdy = PMLy_stretch_PS .* ifft(1i*ky_PS.*fft(vx_PS,[],1),[],1,'symmetric');
            dvydy = PMLy_stretch_PS .* ifft(1i*ky_PS.*fft(vy_PS,[],1),[],1,'symmetric');

            % compute spatial derivatives of physical fluxes
            fluxpa_x = rcc_PS.*dvxdx + dpadx.*ux_PS;
            fluxvx_x = dpadx/rho0 + dvxdx.*ux_PS;
            fluxvy_x = dvydx.*ux_PS;
            fluxpa_y = rcc_PS.*dvydy + dpady.*uy_PS;
            fluxvx_y = dvxdy.*uy_PS;
            fluxvy_y = dpady/rho0 + dvydy.*uy_PS;

            % compute time derivatives
            dpadt = -fluxpa_x -fluxpa_y ;
            dvxdt = -fluxvx_x -fluxvx_y -(vx_PS.*duxdx_PS+vy_PS.*duxdy_PS);
            dvydt = -fluxvy_x -fluxvy_y -(vx_PS.*duydx_PS+vy_PS.*duydy_PS);

            % modify time derivatives to account for PMLs
            dpadt(PMLx_vmap_PS) = dpadt(PMLx_vmap_PS) - PMLx_sigma_PS.*PMLx_pa_PS;
            dvxdt(PMLx_vmap_PS) = dvxdt(PMLx_vmap_PS) - PMLx_sigma_PS.*PMLx_vx_PS;
            dvydt(PMLx_vmap_PS) = dvydt(PMLx_vmap_PS) - PMLx_sigma_PS.*PMLx_vy_PS;
            dpadt(PMLy_vmap_PS) = dpadt(PMLy_vmap_PS) - PMLy_sigma_PS.*PMLy_pa_PS;
            dvxdt(PMLy_vmap_PS) = dvxdt(PMLy_vmap_PS) - PMLy_sigma_PS.*PMLy_vx_PS;
            dvydt(PMLy_vmap_PS) = dvydt(PMLy_vmap_PS) - PMLy_sigma_PS.*PMLy_vy_PS;

            % compute time derivatives of PML auxiliary variables
            PMLx_dpadt = -fluxpa_x(PMLx_vmap_PS)-PMLx_sigma_PS.*PMLx_pa_PS;
            PMLx_dvxdt = -fluxvx_x(PMLx_vmap_PS)-PMLx_sigma_PS.*PMLx_vx_PS;
            PMLx_dvydt = -fluxvy_x(PMLx_vmap_PS)-PMLx_sigma_PS.*PMLx_vy_PS;
            PMLy_dpadt = -fluxpa_y(PMLy_vmap_PS)-PMLy_sigma_PS.*PMLy_pa_PS;
            PMLy_dvxdt = -fluxvx_y(PMLy_vmap_PS)-PMLy_sigma_PS.*PMLy_vx_PS;
            PMLy_dvydt = -fluxvy_y(PMLy_vmap_PS)-PMLy_sigma_PS.*PMLy_vy_PS;

            % space-time transformation to stabilize x-PMLs in the presence of wind
            if PMLx_beta>0
                PMLx_dpadt = PMLx_dpadt+PMLx_beta*(rcc_PS(PMLx_vmap_PS).*dvxdt(PMLx_vmap_PS)+dpadt(PMLx_vmap_PS).*ux_PS(PMLx_vmap_PS));
                PMLx_dvxdt = PMLx_dvxdt+PMLx_beta*(dpadt(PMLx_vmap_PS)/rho0 + dvxdt(PMLx_vmap_PS).*ux_PS(PMLx_vmap_PS));
                PMLx_dvydt = PMLx_dvydt+PMLx_beta*(dvydt(PMLx_vmap_PS).*ux_PS(PMLx_vmap_PS));
            end

            % time integration
            PMLx_pa_PS = PMLx_pa0_PS + dt_PS*RKcoef_PS(rkstep)*PMLx_dpadt;
            PMLx_vx_PS = PMLx_vx0_PS + dt_PS*RKcoef_PS(rkstep)*PMLx_dvxdt;
            PMLx_vy_PS = PMLx_vy0_PS + dt_PS*RKcoef_PS(rkstep)*PMLx_dvydt;

            PMLy_pa_PS = PMLy_pa0_PS + dt_PS*RKcoef_PS(rkstep)*PMLy_dpadt;
            PMLy_vx_PS = PMLy_vx0_PS + dt_PS*RKcoef_PS(rkstep)*PMLy_dvxdt;
            PMLy_vy_PS = PMLy_vy0_PS + dt_PS*RKcoef_PS(rkstep)*PMLy_dvydt;

            pa_PS = pa0_PS + dt_PS * RKcoef_PS(rkstep) * dpadt;
            vx_PS = vx0_PS + dt_PS * RKcoef_PS(rkstep) * dvxdt;
            vy_PS = vy0_PS + dt_PS * RKcoef_PS(rkstep) * dvydt;
            time_PS = time0_PS + dt_PS * RKcoef_PS(rkstep);

        end % PS time integration loop

        % PS filtering
        pa_PS = ifft2(filter_xy.*fft2(pa_PS),'symmetric');
        vx_PS = ifft2(filter_xy.*fft2(vx_PS),'symmetric');
        vy_PS = ifft2(filter_xy.*fft2(vy_PS),'symmetric');

    end

    for DG_cycles = 1:dt_factor % DG time loop
        % n = (n_PS-1)*dt_factor + DG_cycles % DG time iteration
        for rkstep=1:length(rkf84a) % DG time integration loop

            % save external solution across shared intersections
            pap = pa(vmapP);
            vxp = vx(vmapP);
            vyp = vy(vmapP);

            % compute incoming characteristic waves on boundaries
            qcinc = ( pa(vmapB) + rho0*c_eff(vmapB).*( nx(mapB).*vx(vmapB) + ny(mapB).*vy(vmapB)) )/2;
            qcv = ny(mapB).*vx(vmapB) - nx(mapB).*vy(vmapB); % vorticity wave

            % set boundary condition for vorticity wave
            qcv( (nx(mapB).*ux(vmapB)+ny(mapB).*uy(vmapB))<0 ) = 0;

            % compute reflected characteric wave from reflection coefficient
            for ibc=BCcodesInst % only instantenuous response
                qcref(bmapsBC{ibc}) = BCparams{ibc}.Rinf*qcinc(bmapsBC{ibc});
            end

            for ibc=BCcodesReal % real poles (order 1)
                qcref(bmapsBC{ibc}) = BCparams{ibc}.Rinf*qcinc(bmapsBC{ibc}) ...
                    + phi{ibc}*BCparams{ibc}.C;

                % update auxiliary variables
                resphi{ibc} = rkf84a(rkstep)*resphi{ibc} + dt * ...
                    (phi{ibc}.*BCparams{ibc}.lambda.' + qcinc(bmapsBC{ibc}));
                phi{ibc} = phi{ibc} + rkf84b(rkstep)*resphi{ibc};
            end

            for ibc=BCcodesCplx % complex poles (order 2)
                qcref(bmapsBC{ibc}) = BCparams{ibc}.Rinf*qcinc(bmapsBC{ibc}) ...
                    + 2*psir{ibc}*BCparams{ibc}.A + 2*psii{ibc}*BCparams{ibc}.B;

                % update auxiliary variables
                respsir{ibc} = rkf84a(rkstep)*respsir{ibc} + dt * ...
                    (psir{ibc}.*BCparams{ibc}.alpha.' + psii{ibc}.*BCparams{ibc}.beta.' + qcinc(bmapsBC{ibc}));
                respsii{ibc} = rkf84a(rkstep)*respsii{ibc} + dt * ...
                    (psii{ibc}.*BCparams{ibc}.alpha.' - psir{ibc}.*BCparams{ibc}.beta.');

                psir{ibc} = psir{ibc} + rkf84b(rkstep)*respsir{ibc};
                psii{ibc} = psii{ibc} + rkf84b(rkstep)*respsii{ibc};
            end

            % revise external solution on boundaries
            pap(mapB) = qcinc + qcref;
            vxp(mapB) = nx(mapB)./(rho0*c_eff(vmapB)).*(qcinc - qcref) + ny(mapB).*qcv;
            vyp(mapB) = ny(mapB)./(rho0*c_eff(vmapB)).*(qcinc - qcref) - nx(mapB).*qcv;

            % jump differences across shared intersections
            dpa = (pa(vmapM)-pap)/2;
            dvx = (vx(vmapM)-vxp)/2;
            dvy = (vy(vmapM)-vyp)/2;

            % mean across shared intersections
            mpa = (pa(vmapM)+pap)/2;
            mvx = (vx(vmapM)+vxp)/2;
            mvy = (vy(vmapM)+vyp)/2;

            % compute upwind numerical and physical fluxes on intersections
            numfluxpa = c_eff_face.*dpa + rccnx.*mvx + rccny.*mvy + un.*mpa + rcnxun.*dvx + rcnyun.*dvy;
            numfluxvx = nxrho.*mpa + cnxnx.*dvx + cnxny.*dvy + nx.*un.*dpa./(c_eff_face*rho0) + mvx.*un + (dvx -nx.*nx.*dvx -nx.*ny.*dvy).*abs(un);
            numfluxvy = nyrho.*mpa + cnxny.*dvx + cnyny.*dvy + ny.*un.*dpa./(c_eff_face*rho0) + mvy.*un + (dvy -nx.*ny.*dvx -ny.*ny.*dvy).*abs(un);

            physfluxpa = rccnx.*vx(vmapM) + rccny.*vy(vmapM) + un.*pa(vmapM);
            physfluxvx = nxrho.*pa(vmapM) + un.*vx(vmapM);
            physfluxvy = nyrho.*pa(vmapM) + un.*vy(vmapM);

            % compute spatial derivatives
            dpadx = rx.*(Dr*pa)+sx.*(Ds*pa);
            dvxdx = rx.*(Dr*vx)+sx.*(Ds*vx);
            dvydx = rx.*(Dr*vy)+sx.*(Ds*vy);
            dpady = ry.*(Dr*pa)+sy.*(Ds*pa);
            dvxdy = ry.*(Dr*vx)+sy.*(Ds*vx);
            dvydy = ry.*(Dr*vy)+sy.*(Ds*vy);

            % compute spatial derivatives of physical fluxes over the whole domain
            fluxpa_x = rcc.*dvxdx + dpadx.*ux;
            fluxvx_x = dpadx/rho0 + dvxdx.*ux;
            fluxvy_x = dvydx.*ux;
            fluxpa_y = rcc.*dvydy + dpady.*uy;
            fluxvx_y = dvxdy.*uy;
            fluxvy_y = dpady/rho0 + dvydy.*uy;

            % compute time derivatives
            dpadt = -fluxpa_x -fluxpa_y + LIFT*(Fscale.*(physfluxpa-numfluxpa));
            dvxdt = -fluxvx_x -fluxvx_y -(vx.*duxdx+vy.*duxdy) + LIFT*(Fscale.*(physfluxvx-numfluxvx));
            dvydt = -fluxvy_x -fluxvy_y -(vx.*duydx+vy.*duydy) + LIFT*(Fscale.*(physfluxvy-numfluxvy));

            % modify time derivatives to account for PMLs
            dpadt(PMLx_vmapM) = dpadt(PMLx_vmapM) - PMLx_sigma.*PMLx_pa;
            dvxdt(PMLx_vmapM) = dvxdt(PMLx_vmapM) - PMLx_sigma.*PMLx_vx;
            dvydt(PMLx_vmapM) = dvydt(PMLx_vmapM) - PMLx_sigma.*PMLx_vy;
            dpadt(PMLy_vmapM) = dpadt(PMLy_vmapM) - PMLy_sigma.*PMLy_pa;
            dvxdt(PMLy_vmapM) = dvxdt(PMLy_vmapM) - PMLy_sigma.*PMLy_vx;
            dvydt(PMLy_vmapM) = dvydt(PMLy_vmapM) - PMLy_sigma.*PMLy_vy;

            PMLx_dpadt = -fluxpa_x(PMLx_vmapM)-PMLx_sigma.*PMLx_pa;
            PMLx_dvxdt = -fluxvx_x(PMLx_vmapM)-PMLx_sigma.*PMLx_vx;
            PMLx_dvydt = -fluxvy_x(PMLx_vmapM)-PMLx_sigma.*PMLx_vy;
            PMLy_dpadt = -fluxpa_y(PMLy_vmapM)-PMLy_sigma.*PMLy_pa;
            PMLy_dvxdt = -fluxvx_y(PMLy_vmapM)-PMLy_sigma.*PMLy_vx;
            PMLy_dvydt = -fluxvy_y(PMLy_vmapM)-PMLy_sigma.*PMLy_vy;

            % space-time transformation to stabilize x-PMLs in the presence of wind
            if PMLx_beta>0
                PMLx_dpadt = PMLx_dpadt+PMLx_beta*(rcc(PMLx_vmapM).*dvxdt(PMLx_vmapM)+dpadt(PMLx_vmapM).*ux(PMLx_vmapM));
                PMLx_dvxdt = PMLx_dvxdt+PMLx_beta*(dpadt(PMLx_vmapM)/rho0 + dvxdt(PMLx_vmapM).*ux(PMLx_vmapM));
                PMLx_dvydt = PMLx_dvydt+PMLx_beta*(dvydt(PMLx_vmapM).*ux(PMLx_vmapM));
            end

            % increment R-K residuals
            PMLx_respa = rkf84a(rkstep)*PMLx_respa + dt*PMLx_dpadt;
            PMLx_resvx = rkf84a(rkstep)*PMLx_resvx + dt*PMLx_dvxdt;
            PMLx_resvy = rkf84a(rkstep)*PMLx_resvy + dt*PMLx_dvydt;
            PMLy_respa = rkf84a(rkstep)*PMLy_respa + dt*PMLy_dpadt;
            PMLy_resvx = rkf84a(rkstep)*PMLy_resvx + dt*PMLy_dvxdt;
            PMLy_resvy = rkf84a(rkstep)*PMLy_resvy + dt*PMLy_dvydt;

            respa = rkf84a(rkstep)*respa + dt*dpadt;
            resvx = rkf84a(rkstep)*resvx + dt*dvxdt;
            resvy = rkf84a(rkstep)*resvy + dt*dvydt;

            % update fields
            PMLx_pa = PMLx_pa + rkf84b(rkstep)*PMLx_respa;
            PMLx_vx = PMLx_vx + rkf84b(rkstep)*PMLx_resvx;
            PMLx_vy = PMLx_vy + rkf84b(rkstep)*PMLx_resvy;
            PMLy_pa = PMLy_pa + rkf84b(rkstep)*PMLy_respa;
            PMLy_vx = PMLy_vx + rkf84b(rkstep)*PMLy_resvx;
            PMLy_vy = PMLy_vy + rkf84b(rkstep)*PMLy_resvy;

            pa = pa + rkf84b(rkstep)*respa;
            vx = vx + rkf84b(rkstep)*resvx;
            vy = vy + rkf84b(rkstep)*resvy;

        end % DG time integration loop
        time = time + dt; % increment DG time

        % apply modal filter (to try to fix some unstable cases)
        PMLx_pa = Filtmod*PMLx_pa;
        PMLx_vx = Filtmod*PMLx_vx;
        PMLx_vy = Filtmod*PMLx_vy;
        PMLy_pa = Filtmod*PMLy_pa;
        PMLy_vx = Filtmod*PMLy_vx;
        PMLy_vy = Filtmod*PMLy_vy;
        pa = Filtmod*pa;
        vx = Filtmod*vx;
        vy = Filtmod*vy;

    end % DG time loop


    % copy PS to DG
    if ~op_DG_only
        % for corners nodes, no need for interpolation
        pa(map_DG_PS2DG_C) = pa_PS(map_PS_PS2DG_C);
        vx(map_DG_PS2DG_C) = vx_PS(map_PS_PS2DG_C);
        vy(map_DG_PS2DG_C) = vy_PS(map_PS_PS2DG_C);

        % compute spatial FFT of PS fields along x-axis
        pa_fftX = fft(pa_PS,[],2);
        vx_fftX = fft(vx_PS,[],2);
        vy_fftX = fft(vy_PS,[],2);

        % compute spatial FFT of PS fields along y-axis
        pa_fftY = fft(pa_PS,[],1);
        vx_fftY = fft(vx_PS,[],1);
        vy_fftY = fft(vy_PS,[],1);

        for nn = 1:N_dist1
            % spectral interpolation of PS fields along x-axis
            pa_interp_x = ifft(exp(1i*kx_PS*int_dist_vec1(nn)).*pa_fftX,[],2,'symmetric');
            vx_interp_x = ifft(exp(1i*kx_PS*int_dist_vec1(nn)).*vx_fftX,[],2,'symmetric');
            vy_interp_x = ifft(exp(1i*kx_PS*int_dist_vec1(nn)).*vy_fftX,[],2,'symmetric');

            % copy interpolated fields to DG boundary nodes (x-axis)
            pa(map_DG_PS2DG_Bx{nn}) = pa_interp_x(map_PS_PS2DG_Bx{nn});
            vx(map_DG_PS2DG_Bx{nn}) = vx_interp_x(map_PS_PS2DG_Bx{nn});
            vy(map_DG_PS2DG_Bx{nn}) = vy_interp_x(map_PS_PS2DG_Bx{nn});

            % spectral interpolation of PS fields along y-axis
            pa_interp_y = ifft(exp(1i*ky_PS*int_dist_vec1(nn)).*pa_fftY,[],1,'symmetric');
            vx_interp_y = ifft(exp(1i*ky_PS*int_dist_vec1(nn)).*vx_fftY,[],1,'symmetric');
            vy_interp_y = ifft(exp(1i*ky_PS*int_dist_vec1(nn)).*vy_fftY,[],1,'symmetric');

            % copy interpolated fields to DG boundary nodes (y-axis)
            pa(map_DG_PS2DG_By{nn}) = pa_interp_y(map_PS_PS2DG_By{nn});
            vx(map_DG_PS2DG_By{nn}) = vx_interp_y(map_PS_PS2DG_By{nn});
            vy(map_DG_PS2DG_By{nn}) = vy_interp_y(map_PS_PS2DG_By{nn});

            for mm = 1:N_dist1
                % spectral interpolation of PS fields for both axes
                pa_interp_xy = ifft(exp(1i.*ky_PS.*int_dist_vec1(mm)).*fft(pa_interp_x,[],1),[],1,'symmetric');
                vx_interp_xy = ifft(exp(1i.*ky_PS.*int_dist_vec1(mm)).*fft(vx_interp_x,[],1),[],1,'symmetric');
                vy_interp_xy = ifft(exp(1i.*ky_PS.*int_dist_vec1(mm)).*fft(vy_interp_x,[],1),[],1,'symmetric');

                % copy interpolated fields to DG inner nodes
                pa(map_DG_PS2DG_Ixy{mm,nn}) = pa_interp_xy(map_PS_PS2DG_Ixy{mm,nn});
                vx(map_DG_PS2DG_Ixy{mm,nn}) = vx_interp_xy(map_PS_PS2DG_Ixy{mm,nn});
                vy(map_DG_PS2DG_Ixy{mm,nn}) = vy_interp_xy(map_PS_PS2DG_Ixy{mm,nn});

            end

        end

        % copy DG to PS
        pa_PS = gather(pa_PS);
        vx_PS = gather(vx_PS);
        vy_PS = gather(vy_PS);
        pa_PS(map_PS_DG2PS) = gather(pa(map_DG_DG2PS));
        vx_PS(map_PS_DG2PS) = gather(vx(map_DG_DG2PS));
        vy_PS(map_PS_DG2PS) = gather(vy(map_DG_DG2PS));
        pa_PS = gpuArray(pa_PS);
        vx_PS = gpuArray(vx_PS);
        vy_PS = gpuArray(vy_PS);

    end

    % break in case of instability
    simulation_unstable = false;
    if max(abs(pa(:)))>10*pa_max_ini
        simulation_unstable = true;
    end
    if ~op_DG_only
        if max(abs(pa_PS(:)))>10*pa_max_ini
            simulation_unstable = true;
        end
    end
    if simulation_unstable
        beep, pause(0.01), beep
        warning('simulation unstable!')
        break
    end

    % store pressure signals at recording locations
    if op_record_DG
        for irec=1:length(rec_id)
            pa_rec(irec,n_PS) = rec_interp(irec,:)*pa(:,rec_id(irec));
        end
    end

    if op_record_PS
        for irec=1:length(rec_idx_PS)
            pa_rec_PS(irec,n_PS) = pa_PS(rec_idy_PS(irec),rec_idx_PS(irec));
        end
    end

    % snapshot of pressure field over the whole DG domain
    if mod(n_PS,op_snapshot)==0
        if op_plot

            try
                if ~eq(gcf,fig_plt)
                    figure(fig_plt)
                end
            catch
                fig_plt = figure(667);
            end

            plottingScript;
            colorbar
            drawnow
            if op_print_snapshots
                op_print_snapshotf(n_PS,output_name);
            end

        end % op_plot

        % write spatial pressure field to disk
        if op_record_space
            filename=sprintf([output_name '_%d.mat'],n_PS);
            disp(['Writing ' filename ' to disk'])
            pressure=gather(pa);
            save(fullfile('outputs',filename),'pressure');
        end

    end % snapshot

    % save recorded pressure signals (this is done at each iteration, instead of
    % at the end of the simulation, in case something goes wrong)
    if op_record_DG
        parec=gather(pa_rec);
        dtrec=gather(dt_PS);
        xrec=gather(recx); yrec=gather(recy);
        save(fullfile('outputs',[output_name '_rec.mat']),'parec','xrec','yrec','dtrec');
    end

    if op_record_PS
        parec=gather(pa_rec_PS);
        dtrec=gather(dt_PS);
        xrec=gather(recx_PS); yrec=gather(recy_PS);
        save(fullfile('outputs',[output_name '_rec_PS.mat']),'parec','xrec','yrec','dtrec');
    end

end % PS time loop
toc

% elapsed physical time
time=gather(time);
