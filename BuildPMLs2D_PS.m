% Note: initial spatial source distribution must already be defined to initialize the PMLs

% Compute position of internal PML boundaries
PMLxb_ext_PS = min(x_PS(:));
PMLxt_ext_PS = max(x_PS(:));
PMLxb_int_PS = PMLxb_ext_PS+PMLxb_thickness_PS;
PMLxt_int_PS = PMLxt_ext_PS-PMLxt_thickness_PS;

PMLyb_ext_PS = min(y_PS(:));
PMLyt_ext_PS = max(y_PS(:));
PMLyb_int_PS = PMLyb_ext_PS+PMLyb_thickness_PS;
PMLyt_int_PS = PMLyt_ext_PS-PMLyt_thickness_PS;

if (PMLxt_int_PS<PMLxb_int_PS) || (PMLyt_int_PS<PMLyb_int_PS)
    error('PMLs should not overlap')
end

% find number of PML points along each directions
PMLxb_Nx_PS = length(find(xvec_PS<PMLxb_int_PS-NODETOL));
PMLxt_Nx_PS = length(find(xvec_PS>PMLxt_int_PS+NODETOL));
PMLxb_Ny_PS = length(yvec_PS);
PMLxt_Ny_PS = length(yvec_PS);
PMLyb_Ny_PS = length(find(yvec_PS<PMLyb_int_PS-NODETOL));
PMLyt_Ny_PS = length(find(yvec_PS>PMLyt_int_PS+NODETOL));
PMLyb_Nx_PS = length(xvec_PS);
PMLyt_Nx_PS = length(xvec_PS);

% find points located in the PMLs
PMLxb_vmap_PS = find(x_PS<PMLxb_int_PS-NODETOL);
% $$$ if any(PMLxb_vmap_PS) % add 1 element to simplify hybrid maps
% $$$     PMLxb_vmap_PS = find(x_PS<PMLxb_int_PS+dx_PS-NODETOL);
% $$$ end
PMLxt_vmap_PS = find(x_PS>PMLxt_int_PS+NODETOL);
% $$$ if any(PMLxt_vmap_PS)
% $$$     PMLxt_vmap_PS = find(x_PS>PMLxt_int_PS-dx_PS+NODETOL);
% $$$ end
PMLyb_vmap_PS = find(y_PS<PMLyb_int_PS-NODETOL);
PMLyt_vmap_PS = find(y_PS>PMLyt_int_PS+NODETOL);

% $$$ % in case hybrid model
% $$$ if any(PMLxb_vmap_PS) % add 1 element for hybrid maps
% $$$     PMLxb_Nx_PS = length(find(xvec_PS<PMLxb_int_PS+dx_PS-NODETOL));
% $$$ end
% $$$ if any(PMLxt_vmap_PS)
% $$$     PMLxt_Nx_PS = length(find(xvec_PS>PMLxt_int_PS-dx_PS+NODETOL));
% $$$ end

% wavenumber PML vectors (e.g., for filtering purposes)
PMLxb_kxvec_PS = wavenumber_vector(PMLxb_Nx_PS, dx_PS);
PMLxb_kyvec_PS = wavenumber_vector(PMLxb_Ny_PS, dy_PS);
PMLxt_kxvec_PS = wavenumber_vector(PMLxt_Nx_PS, dx_PS);
PMLxt_kyvec_PS = wavenumber_vector(PMLxt_Ny_PS, dy_PS);
PMLyb_kxvec_PS = wavenumber_vector(PMLyb_Nx_PS, dx_PS);
PMLyb_kyvec_PS = wavenumber_vector(PMLyb_Ny_PS, dy_PS);
PMLyt_kxvec_PS = wavenumber_vector(PMLyt_Nx_PS, dx_PS);
PMLyt_kyvec_PS = wavenumber_vector(PMLyt_Ny_PS, dy_PS);

% corresponding matrices
PMLxb_kx_PS = repmat(PMLxb_kxvec_PS,PMLxb_Ny_PS,1);
PMLxb_ky_PS = repmat(PMLxb_kyvec_PS,PMLxb_Nx_PS,1).';
PMLxt_kx_PS = repmat(PMLxt_kxvec_PS,PMLxt_Ny_PS,1);
PMLxt_ky_PS = repmat(PMLxt_kyvec_PS,PMLxt_Nx_PS,1).';
PMLyb_kx_PS = repmat(PMLyb_kxvec_PS,PMLyb_Ny_PS,1);
PMLyb_ky_PS = repmat(PMLyb_kyvec_PS,PMLyb_Nx_PS,1).';
PMLyt_kx_PS = repmat(PMLyt_kxvec_PS,PMLyt_Ny_PS,1);
PMLyt_ky_PS = repmat(PMLyt_kyvec_PS,PMLyt_Nx_PS,1).';

% compute relative distance to external boundary
PMLxb_rdist_PS = PMLxb_int_PS - x_PS(PMLxb_vmap_PS);
PMLxt_rdist_PS = x_PS(PMLxt_vmap_PS) - PMLxt_int_PS;
PMLyb_rdist_PS = PMLyb_int_PS - y_PS(PMLyb_vmap_PS);
PMLyt_rdist_PS = y_PS(PMLyt_vmap_PS) - PMLyt_int_PS;
PMLxb_rdist_PS = max(PMLxb_rdist_PS,0);
PMLxt_rdist_PS = max(PMLxt_rdist_PS,0);
PMLyb_rdist_PS = max(PMLyb_rdist_PS,0);
PMLyt_rdist_PS = max(PMLyt_rdist_PS,0);

% compute absorption profiles
PMLxb_sigma_PS = profile_sigma(PMLxb_rdist_PS, PMLxb_sigma_max_PS, PMLxb_pwr_PS, PMLxb_thickness_PS);
PMLxt_sigma_PS = profile_sigma(PMLxt_rdist_PS, PMLxt_sigma_max_PS, PMLxt_pwr_PS, PMLxt_thickness_PS);
PMLyb_sigma_PS = profile_sigma(PMLyb_rdist_PS, PMLyb_sigma_max_PS, PMLyb_pwr_PS, PMLyb_thickness_PS);
PMLyt_sigma_PS = profile_sigma(PMLyt_rdist_PS, PMLyt_sigma_max_PS, PMLyt_pwr_PS, PMLyt_thickness_PS);

% compute stretch profiles
PMLxb_stretch_PS = profile_stretch(PMLxb_rdist_PS, PMLxb_stretch_max_PS, PMLxb_pwr_PS, PMLxb_thickness_PS);
PMLxt_stretch_PS = profile_stretch(PMLxt_rdist_PS, PMLxt_stretch_max_PS, PMLxt_pwr_PS, PMLxt_thickness_PS);
PMLyb_stretch_PS = profile_stretch(PMLyb_rdist_PS, PMLyb_stretch_max_PS, PMLyb_pwr_PS, PMLyb_thickness_PS);
PMLyt_stretch_PS = profile_stretch(PMLyt_rdist_PS, PMLyt_stretch_max_PS, PMLyt_pwr_PS, PMLyt_thickness_PS);

% merge the PML vectors
PMLx_vmap_PS = [PMLxb_vmap_PS; PMLxt_vmap_PS];
PMLy_vmap_PS = [PMLyb_vmap_PS; PMLyt_vmap_PS];

PMLx_sigma_PS = [PMLxb_sigma_PS.*PMLxb_stretch_PS; PMLxt_sigma_PS.*PMLxt_stretch_PS];
PMLy_sigma_PS = [PMLyb_sigma_PS.*PMLyb_stretch_PS; PMLyt_sigma_PS.*PMLyt_stretch_PS];
% $$$ PMLx_sigma_PS = reshape(PMLx_sigma_PS, [length(yvec_PS),PMLxb_Nx_PS+PMLxt_Nx_PS]);
% $$$ PMLy_sigma_PS = reshape(PMLy_sigma_PS, [PMLyb_Ny_PS+PMLyt_Ny_PS],length(xvec_PS));

PMLx_stretch_PS = ones(size(x_PS));
PMLx_stretch_PS(PMLx_vmap_PS) = [PMLxb_stretch_PS; PMLxt_stretch_PS];
PMLy_stretch_PS = ones(size(y_PS));
PMLy_stretch_PS(PMLy_vmap_PS) = [PMLyb_stretch_PS; PMLyt_stretch_PS];

% build local maps so that PMLx_vmap(PMLxb_lvmap)==PMLxb_vmap (e.g., for spectral filtering)
[~,~,PMLxb_map_PS] = intersect(PMLxb_vmap_PS,PMLx_vmap_PS,'stable');
[~,~,PMLxt_map_PS] = intersect(PMLxt_vmap_PS,PMLx_vmap_PS,'stable');
[~,~,PMLyb_map_PS] = intersect(PMLyb_vmap_PS,PMLy_vmap_PS,'stable');
[~,~,PMLyt_map_PS] = intersect(PMLyt_vmap_PS,PMLy_vmap_PS,'stable');

PMLxb_map_PS = reshape(PMLxb_map_PS, [PMLxb_Ny_PS PMLxb_Nx_PS]);
PMLxt_map_PS = reshape(PMLxt_map_PS, [PMLxt_Ny_PS PMLxt_Nx_PS]);
PMLyb_map_PS = reshape(PMLyb_map_PS, [PMLyb_Ny_PS PMLyb_Nx_PS]);
PMLyt_map_PS = reshape(PMLyt_map_PS, [PMLyt_Ny_PS PMLyt_Nx_PS]);


%%
% initialize PML auxiliary variables (not entirely finished, but will work
% fine)

PMLx_pa_PS = -pa_PS(PMLx_vmap_PS).*PMLx_sigma_PS/1;
PMLx_vx_PS = -vx_PS(PMLx_vmap_PS).*PMLx_sigma_PS/1;
PMLx_vy_PS = -vy_PS(PMLx_vmap_PS).*PMLx_sigma_PS/1;

PMLy_pa_PS = -pa_PS(PMLy_vmap_PS).*PMLy_sigma_PS/1;
PMLy_vx_PS = -vx_PS(PMLy_vmap_PS).*PMLy_sigma_PS/1;
PMLy_vy_PS = -vy_PS(PMLy_vmap_PS).*PMLy_sigma_PS/1;

PMLx_pa0_PS = zeros(size(PMLx_vmap_PS));
PMLx_vx0_PS = zeros(size(PMLx_vmap_PS));
PMLx_vy0_PS = zeros(size(PMLx_vmap_PS));

PMLy_pa0_PS = zeros(size(PMLy_vmap_PS));
PMLy_vy0_PS = zeros(size(PMLy_vmap_PS));
PMLy_vy0_PS = zeros(size(PMLy_vmap_PS));

% clear variables not useful anymore
% $$$ clear PMLxb_vmap_PS PMLxt_vmap_PS
% $$$ clear PMLyb_vmap_PS PMLyt_vmap_PS
% $$$ clear PMLxb_rdist_PS PMLxt_rdist_PS
% $$$ clear PMLyb_rdist_PS PMLyt_rdist_PS
% $$$ clear PMLxb_sigma_PS PMLxt_sigma_PS
% $$$ clear PMLyb_sigma_PS PMLyt_sigma_PS
% $$$ clear PMLxb_stretch_PS PMLxt_stretch_PS
% $$$ clear PMLyb_stretch_PS PMLyt_stretch_PS
