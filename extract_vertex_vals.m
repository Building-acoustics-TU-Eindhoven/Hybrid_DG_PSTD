% disable GPU computations
try
    VX = gather(VX);
    VY = gather(VY);
    x = gather(x);
    y = gather(y);
    pa = gather(pa);
end

% convert to linear indexing
xx = x(:);
yy = y(:);
pp = pa(:);

% new acoustic vector
pa_v = zeros(size(VX));

% calculate all pairwise distances
tic
changes = zeros(Np*K,Nv);
for ii = 1:Nv
    changes(:,ii) = (xx-VX(ii)).^2 + (yy-VY(ii)).^2;
    % the square (Euclidean) distance was used, to save on sqrt calculations.
end

% find distances less than a given tolerance
[jj,kk] = find(changes < NODETOL);

for ii=1:Nv

    % find positions corresponding to the same vertex
    idx = find(kk==ii);

    % copy acoustic values
    pa_v(ii) = mean(pp(jj(idx)));

end
toc


%% plot (only works for structured quadrangles)

% round vertices to mesh tolerance for easier processing
VXtmp = round(VX/NODETOL)*NODETOL;
VYtmp = round(VY/NODETOL)*NODETOL;

% grid layout to be used for plotting
xxx = unique(VXtmp);
yyy = unique(VYtmp);
[XXX,YYY]=meshgrid(xxx,yyy);

% Now, we want to reshape VX, VY, and pa_v to a similar layout
VX2 = zeros(size(XXX));
VY2 = zeros(size(XXX));
pa_v2 = zeros(size(XXX));

for ii=1:size(XXX,1)
    for jj=1:size(XXX,2)

        % find matching vertices
        idx = find(VXtmp==XXX(ii,jj) & VYtmp==YYY(ii,jj));

        % copy data
        VX2(ii,jj) = VX(idx);
        VY2(ii,jj) = VY(idx);
        pa_v2(ii,jj) = pa_v(idx);

    end
end

figure
pcolor(VX2,VY2,pa_v2)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
axis image
colorbar