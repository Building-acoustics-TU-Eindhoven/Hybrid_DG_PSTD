function [V2D] = Vandermonde2D_quad(N, r, s);

% function [V2D] = Vandermonde2Dquad(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix,  V_{ij} = phi_j(r_i)phi_j(s_i);

V2D = zeros(length(r(:)),(N+1)*(N+1));

% build the Vandermonde matrix
sk = 1;
for i=0:N
  for j=0:N
    h1 = JacobiP(r(:),0,0,i);
    h2 = JacobiP(s(:),0,0,j);
    V2D(:,sk) = h1.*h2;
    sk = sk+1;
  end
end
return
