function filter = spectral_filter_2D(kx, ky, ppw, alpha, beta, kmax, kc)
% SPECTRAL_FILTER_2D: build a spectral filter to remove the short wave
% components along the x and y axes of an acoustic field (based on
% 10.1016/j.jcp.2017.07.046). The filter can e.g. be applied with
%
% pa = ifft2(filter.*fft2(pa),'symmetric');
%
% kx & ky: wavenumbers along the x and y axes
% ppw, alpha, beta: filter parameters, see paper
% kmax: Nyquist–Shannon wavenumber (e.g., pi/dx)
% kc: cut-off wavenumber

cond1 = find( abs(kx)>kc );
cond2 = find( abs(kx)>kc & abs(ky)<=kc );
cond3 = find( abs(kx)<=kc & abs(ky)>kc );
cond4 = find( abs(kx)<=kc & abs(ky)<=kc );

filter = ones(size(kx));
filter(cond1) = exp( -alpha*log(10).*( (abs(kx(cond1))-kc)./(kmax-kc) ).^(2*beta) + ...
                     -alpha*log(10).*( (abs(ky(cond1))-kc)./(kmax-kc) ).^(2*beta) ...
                        );
filter(cond2) = exp( -alpha*log(10).*( (abs(kx(cond2))-kc)./(kmax-kc) ).^(2*beta) );
filter(cond3) = exp( -alpha*log(10).*( (abs(ky(cond3))-kc)./(kmax-kc) ).^(2*beta) );
filter(cond4) = 1;
