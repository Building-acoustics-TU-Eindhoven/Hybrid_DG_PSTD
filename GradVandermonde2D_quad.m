function [V2Dr,V2Ds] = GradVandermonde2D_quad(N,r,s)

% function [V2Dr,V2Ds] = GradVandermonde2Dquad(N,r,s)
% Purpose : Initialize the gradient of the modal basis (i,j) at (r,s) at order N

V2Dr = zeros(length(r),(N+1)*(N+1));
V2Ds = zeros(length(r),(N+1)*(N+1));

% Initialize matrices
sk = 1;
for i=0:N
  for j=0:N
    [V2Dr(:,sk),V2Ds(:,sk)] = GradSimplex2DP_quad(r,s,i,j);
    sk = sk+1;
  end
end
return;
