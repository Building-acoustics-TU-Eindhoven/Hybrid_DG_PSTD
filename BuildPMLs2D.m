% Note: initial spatial source distribution must already be defined to initialize the PMLs

PMLxb_stretch_max = 1;
PMLxt_stretch_max = 1;
PMLyb_stretch_max = 1;
PMLyt_stretch_max = 1;

%profile_stretch = @(pos,amp,pwr,thick) 1./(1 + (amp-1)*(pos/thick).^pwr );

% Compute position of internal PML boundaries
PMLxb_ext = min(x(:));
PMLxt_ext = max(x(:));
PMLxb_int = PMLxb_ext+PMLxb_thickness;
PMLxt_int = PMLxt_ext-PMLxt_thickness;

PMLyb_ext = min(y(:));
PMLyt_ext = max(y(:));
PMLyb_int = PMLyb_ext+PMLyb_thickness;
PMLyt_int = PMLyt_ext-PMLyt_thickness;

if (PMLxt_int<PMLxb_int) || (PMLyt_int<PMLyb_int)
    error('PMLs should not overlap')
end

% find elements located in the PMLs
[~, PMLxb_emap] = find(x<PMLxb_int-NODETOL);
[~, PMLxt_emap] = find(x>PMLxt_int+NODETOL);
[~, PMLyb_emap] = find(y<PMLyb_int-NODETOL);
[~, PMLyt_emap] = find(y>PMLyt_int+NODETOL);

% if one node of an element is in the PMLs, the whole element is considered
% to be in the PMLs
PMLxb_emap = unique(PMLxb_emap).';
PMLxt_emap = unique(PMLxt_emap).';
PMLyb_emap = unique(PMLyb_emap).';
PMLyt_emap = unique(PMLyt_emap).';

% compute number of PML elements
PMLxb_K = length(PMLxb_emap);
PMLxt_K = length(PMLxt_emap);
PMLyb_K = length(PMLyb_emap);
PMLyt_K = length(PMLyt_emap);
PMLx_K = PMLxb_K + PMLxt_K;
PMLy_K = PMLyb_K + PMLyt_K;

% build PML maps
PMLxb_vmapM = sub2ind(size(x), repmat(1:size(x,1),size(PMLxb_emap)), repelem(PMLxb_emap,size(x,1))).';
PMLxt_vmapM = sub2ind(size(x), repmat(1:size(x,1),size(PMLxt_emap)), repelem(PMLxt_emap,size(x,1))).';
PMLyb_vmapM = sub2ind(size(y), repmat(1:size(y,1),size(PMLyb_emap)), repelem(PMLyb_emap,size(y,1))).';
PMLyt_vmapM = sub2ind(size(y), repmat(1:size(y,1),size(PMLyt_emap)), repelem(PMLyt_emap,size(y,1))).';
PMLxb_vmapM = reshape(PMLxb_vmapM,[Np PMLxb_K]);
PMLxt_vmapM = reshape(PMLxt_vmapM,[Np PMLxt_K]);
PMLyb_vmapM = reshape(PMLyb_vmapM,[Np PMLyb_K]);
PMLyt_vmapM = reshape(PMLyt_vmapM,[Np PMLyt_K]);

% merge the PML vectors
PMLx_emap = [PMLxb_emap.'; PMLxt_emap.'];
PMLy_emap = [PMLyb_emap.'; PMLyt_emap.'];
PMLx_vmapM = [PMLxb_vmapM PMLxt_vmapM];
PMLy_vmapM = [PMLyb_vmapM PMLyt_vmapM];

% construct PML maps
% PML{x,y,z}_vmapM allows to find the PML boundary nodes for matrices
% defined over the whole domain (Omega -> Omega_PML)
% PML{x,y,z}_mapM allows to find the PML boundary nodes for matrices
% defined over the element boundaries (dOmega -> dOmega_PML)
PMLx_mapM = mapM(:,PMLx_emap);
PMLy_mapM = mapM(:,PMLy_emap);
PMLx_mapP = mapP(:,PMLx_emap);
PMLy_mapP = mapP(:,PMLy_emap);

PMLxb_mapM = mapM(:,PMLxb_emap);
PMLxt_mapM = mapM(:,PMLxt_emap);
PMLyb_mapM = mapM(:,PMLyb_emap);
PMLyt_mapM = mapM(:,PMLyt_emap);
% PMLx_mapP = mapP(:,PMLx_emap);
% PMLy_mapP = mapP(:,PMLy_emap);

PMLx_mapB = find(PMLx_mapM==PMLx_mapP);
PMLy_mapB = find(PMLy_mapM==PMLy_mapP);

% compute relative distance to external boundary
PMLxb_rdist = PMLxb_int - x(PMLxb_vmapM);
PMLxt_rdist = x(PMLxt_vmapM) - PMLxt_int;
PMLxb_rdist = max(PMLxb_rdist,0);
PMLxt_rdist = max(PMLxt_rdist,0);
PMLyb_rdist = PMLyb_int - y(PMLyb_vmapM);
PMLyt_rdist = y(PMLyt_vmapM) - PMLyt_int;
PMLyb_rdist = max(PMLyb_rdist,0);
PMLyt_rdist = max(PMLyt_rdist,0);

% % do the same for the element boudnaries \Gamma_e
PMLxb_rdist_face = PMLxb_int - x(PMLxb_mapM);
PMLxt_rdist_face = x(PMLxt_mapM) - PMLxt_int;
PMLxb_rdist_face = max(PMLxb_rdist_face,0);
PMLxt_rdist_face = max(PMLxt_rdist_face,0);
PMLyb_rdist_face = PMLyb_int - y(PMLyb_mapM);
PMLyt_rdist_face = y(PMLyt_mapM) - PMLyt_int;
PMLyb_rdist_face = max(PMLyb_rdist_face,0);
PMLyt_rdist_face = max(PMLyt_rdist_face,0);

% compute absorption profiles
PMLxb_sigma = profile_sigma(PMLxb_rdist, PMLxb_sigma_max, PMLxb_pwr, PMLxb_thickness);
PMLxt_sigma = profile_sigma(PMLxt_rdist, PMLxt_sigma_max, PMLxt_pwr, PMLxt_thickness);
PMLyb_sigma = profile_sigma(PMLyb_rdist, PMLyb_sigma_max, PMLyb_pwr, PMLyb_thickness);
PMLyt_sigma = profile_sigma(PMLyt_rdist, PMLyt_sigma_max, PMLyt_pwr, PMLyt_thickness);

PMLxb_sigma_face = profile_sigma(PMLxb_rdist_face, PMLxb_sigma_max, PMLxb_pwr, PMLxb_thickness);
PMLxt_sigma_face = profile_sigma(PMLxt_rdist_face, PMLxt_sigma_max, PMLxt_pwr, PMLxt_thickness);
PMLyb_sigma_face = profile_sigma(PMLyb_rdist_face, PMLyb_sigma_max, PMLyb_pwr, PMLyb_thickness);
PMLyt_sigma_face = profile_sigma(PMLyt_rdist_face, PMLyt_sigma_max, PMLyt_pwr, PMLyt_thickness);

% PMLxb_stretch = profile_stretch(PMLxb_rdist, PMLxb_stretch_max, PMLxb_pwr, PMLxb_thickness);
% PMLxt_stretch = profile_stretch(PMLxt_rdist, PMLxt_stretch_max, PMLxt_pwr, PMLxt_thickness);
% PMLyb_stretch = profile_stretch(PMLyb_rdist, PMLyb_stretch_max, PMLyb_pwr, PMLyb_thickness);
% PMLyt_stretch = profile_stretch(PMLyt_rdist, PMLyt_stretch_max, PMLyt_pwr, PMLyt_thickness);

PMLxb_stretch_face = profile_stretch(PMLxb_rdist_face, PMLxb_stretch_max, PMLxb_pwr, PMLxb_thickness);
PMLxt_stretch_face = profile_stretch(PMLxt_rdist_face, PMLxt_stretch_max, PMLxt_pwr, PMLxt_thickness);
PMLyb_stretch_face = profile_stretch(PMLyb_rdist_face, PMLyb_stretch_max, PMLyb_pwr, PMLyb_thickness);
PMLyt_stretch_face = profile_stretch(PMLyt_rdist_face, PMLyt_stretch_max, PMLyt_pwr, PMLyt_thickness);


% merge absorption profiles
PMLx_sigma = [PMLxb_sigma PMLxt_sigma];
PMLy_sigma = [PMLyb_sigma PMLyt_sigma];
PMLx_sigma_face = [PMLxb_sigma_face PMLxt_sigma_face];
PMLy_sigma_face = [PMLyb_sigma_face PMLyt_sigma_face];
%PMLx_stretch = [PMLxb_stretch PMLxt_stretch];
%PMLy_stretch = [PMLyb_stretch PMLyt_stretch];
PMLx_stretch_face = [PMLxb_stretch_face PMLxt_stretch_face];
PMLy_stretch_face = [PMLyb_stretch_face PMLyt_stretch_face];


% check that the mesh is structured inside the PMLs
structured_PMLx = all(abs(nx(PMLx_mapM(:)))>1-NODETOL | abs(ny(PMLx_mapM(:)))>1-NODETOL);
structured_PMLy = all(abs(nx(PMLy_mapM(:)))>1-NODETOL | abs(ny(PMLy_mapM(:)))>1-NODETOL);
if ~structured_PMLx || ~structured_PMLy
    warning('The mesh needs structured quadrangles inside the PMLs!')
end


% initialize PML auxiliary variables
PMLx_pa = -pa(PMLx_vmapM).*ones(size(PMLx_vmapM));
PMLx_vx = -vx(PMLx_vmapM).*ones(size(PMLx_vmapM));
PMLx_vy = -vy(PMLx_vmapM).*ones(size(PMLx_vmapM));

PMLy_pa = -pa(PMLy_vmapM).*ones(size(PMLy_vmapM));
PMLy_vx = -vx(PMLy_vmapM).*ones(size(PMLy_vmapM));
PMLy_vy = -vy(PMLy_vmapM).*ones(size(PMLy_vmapM));

PMLx_respa = zeros(size(PMLx_vmapM));
PMLx_resvx = zeros(size(PMLx_vmapM));
PMLx_resvy = zeros(size(PMLx_vmapM));

PMLy_respa = zeros(size(PMLy_vmapM));
PMLy_resvx = zeros(size(PMLy_vmapM));
PMLy_resvy = zeros(size(PMLy_vmapM));

% clear variables not useful anymore
clear PMLxb_emap PMLxt_emap
clear PMLyb_emap PMLyt_emap
clear PMLxb_vmapM PMLxt_vmapM
clear PMLyb_vmapM PMLyt_vmapM
clear PMLxb_rdist PMLxt_rdist
clear PMLyb_rdist PMLyt_rdist
clear PMLxb_sigma PMLxt_sigma
clear PMLyb_sigma PMLyt_sigma
%clear PMLxb_stretch PMLxt_stretch
%clear PMLyb_stretch PMLyt_stretch
%clear PMLx_stretch PMLy_stretch %not used for now
