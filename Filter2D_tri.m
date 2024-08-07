function [F] = Filter2D_quad(Nc,sp,alpha_opt)
% function [F] = Filter2D(Norder,Nc,sp)
% Purpose : Initialize 2D filter matrix of order sp and cutoff Nc

Globals2D;

filterdiag = ones((Norder+1)*(Norder+2)/2,1);

if nargin>=3
    alpha = alpha_opt;
else
    alpha = -log(eps);
    alpha = 8;
end

% build exponential filter
sk = 1;
for i=0:Norder
    for j=0:Norder-i
        if (i+j>=Nc)
            filterdiag(sk) = exp(-alpha*((i+j - Nc)/(Norder-Nc))^sp);
        end
        sk = sk+1;
    end
end

F = Vanderm*diag(filterdiag)*invVanderm;

return;
