function [V2D] = Vandermonde2D_tri(N, r, s);

% function [V2D] = Vandermonde2D(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix,  V_{ij} = phi_j(r_i, s_i);

V2D = zeros(length(r),(N+1)*(N+2)/2);

% Transfer to (a,b) coordinates
[a, b] = rstoab_tri(r, s);

% build the Vandermonde matrix
sk = 1;
for ii=0:N
  for jj=0:N - ii
      V2D(:,sk) = Simplex2DP_tri(a,b,ii,jj);
    sk = sk+1;
  end
end
return;
