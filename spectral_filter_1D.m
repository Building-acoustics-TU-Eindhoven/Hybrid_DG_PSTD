function filter = spectral_filter_1D(kx, ppw, alpha, beta, kmax, kc)
% SPECTRAL_FILTER_1D: build a spectral filter to remove the short wave
% components along one specific axis of an acoustic field (based on
% 10.1016/j.jcp.2017.07.046). The filter can e.g. be applied with
%
% pa = ifft(filter.*fft(pa),'symmetric');
%
% kx: wavenumber along considered axis
% ppw, alpha, beta: filter parameters, see paper
% kmax: Nyquist–Shannon wavenumber (e.g., pi/dx)
% kc: cut-off wavenumber

cond1 = find( abs(kx)>kc );
cond2 = find( abs(kx)<=kc );

filter = ones(size(kx));
filter(cond1) = exp( -alpha*log(10).*( (abs(kx(cond1))-kc)./(kmax-kc) ).^(2*beta) );
filter(cond2) = 1;
