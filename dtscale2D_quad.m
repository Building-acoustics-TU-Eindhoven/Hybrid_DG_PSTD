function dtscale = dtscale2D_quad;

% function dtscale = dtscale2Dquad;
% Purpose : Compute inscribed circle diameter as characteristic
%           for grid to choose timestep
% By Allan P. Engsig-Karup.

Globals2D;

% Find vertex nodes
vx = VX(EToV); vy = VY(EToV);

% Compute estimate of area of quadrilaterals based on two side vectors
if length(vx(:))==4
    dx1 = vx(2) - vx(1);
    dx2 = vx(4) - vx(1);
    dy1 = vy(2) - vy(1);
    dy2 = vy(4) - vy(1);
else
    dx1 = vx(:,2) - vx(:,1);
    dx2 = vx(:,4) - vx(:,1);
    dy1 = vy(:,2) - vy(:,1);
    dy2 = vy(:,4) - vy(:,1);
end
Area = cross([dx1 dy1 dx1*0],[dx2 dy2 dy2*0],2);
Area = max(Area,[],2);
if min(Area)<0
    disp('Some quadrilaterals are reordered anti-clockwise.')
    wer
end

% compute estimated perimeter
P = 2*(sqrt(dx1.^2+dy1.^2) + sqrt(dx2.^2+dy2.^2)); % parallellogram assumed

% compute estimate of scale using ratio between area and perimeter of
% quadrailateral
dtscale = 4*Area./P;
return;
