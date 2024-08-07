function [t, p] = tsearchn_mod(x,elem,xi)
%TSEARCHN_MOD N-D closest simplex search.
% This is a modification of the Matlab function tsearchn to allow the use of
% quadrangles (as well as triangles). Note that this function can be less
% efficient than the original function for large meshes.

%   Copyright 1984-2016 The MathWorks, Inc. + Matthias Cosnefroy 2020

[npt, ndim] = size(xi);                 % Number of points
myeps = builtin('_tsearchntol', ndim);  % Tolerance
myeps = -myeps;
nelem = size(elem,1);                     % Number of simplexes
nnodes = size(elem,2); % Number of geometrical nodes

if size(x,2) ~= ndim
    error(message('MATLAB:tsearchn:InvalidDimensions'))
end

algorithm_type = 2;

t = nan(npt,1);            % Simplex containing corresponding input point
p = nan(npt,nnodes);       % Barycentric coordinates for corresponding input point
X = [ones(size(x,1),1) x]; % Append 1s to vertex matrix
b = [ones(npt,1) xi];      % Append 1s to point matrix

for i = 1:nelem             % Return the largest simplex index for each element
    q = b / X(elem(i,:),:);   % Compute barycentric coordinate of each point
    I = all(q > myeps,2);    % Find simplex where all coordinates are positive
    t(I) = i;                % Set simplex
    p(I,:) = q(I,:);         % Set barycentric coordinates
end
