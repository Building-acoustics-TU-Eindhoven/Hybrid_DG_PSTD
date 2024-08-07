function [V2Dr,V2Ds] = GradVandermonde2D_tri(N,r,s)

% function [V2Dr,V2Ds] = GradVandermonde2D(N,r,s)
% Purpose : Initialize the gradient of the modal basis (i,j) at (r,s) at order N

V2Dr = zeros(length(r),(N+1)*(N+2)/2);
V2Ds = zeros(length(r),(N+1)*(N+2)/2);

% find tensor-product coordinates
[a,b] = rstoab_tri(r,s);

% Initialize matrices
sk = 1;
for ii=0:N
  for jj=0:N-ii
    [V2Dr(:,sk),V2Ds(:,sk)] = GradSimplex2DP_tri(a,b,ii,jj);
    sk = sk+1;
  end
end
return;
