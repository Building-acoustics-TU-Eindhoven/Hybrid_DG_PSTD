% Purpose: declare global variables

% DG parameters
global Np Nfp Norder Kel
global Nfaces EToE EToF EToV
global Vanderm invVanderm
global NODETOL VX VY

NODETOL = 1e-8; % tolerance for mesh node positions

% numerical coefficients
rkf84a = [0 -0.553443129450157 0.0106598757020349 -0.551581288893200 -1.88579037755874 -5.70129574279326 2.11390396566479 -0.533957882667528];
rkf84b = [0.0803793688273695 0.538849745856984 0.0197497440903196 0.0991184129733997 0.746692041106412 1.67958424561889 0.243372806700819 0.142273045900137];
rkf84c = [0 0.0803793688273695 0.321006425033843 0.340850182660466 0.385036482428547 0.504005247753410 0.657897756116854 0.948408762334848];
rkf84g = [1 1/2 1/6 1/24 8.02921837189987e-3 1.10873426499598e-3 9.46273413180222e-5 3.68184991253961e-6];

% Toulorge thesis. Error magnitude: (1) Emag=1dB; (2) Emag=0.1dB; (3) Emag=0.01dB; (4) Emag=0.001dB; (5) Emag=0.0001dB;
CFL_DG_values = [1.19 0.619 0.394 0.276 0.206 0.16 0.127 0.104 0.0873 0.0742;...
                 1.19 0.619 0.394 0.276 0.206 0.16 0.127 0.104 0.0873 0.0742;...
                 1.19 0.619 0.394 0.276 0.206 0.16 0.127 0.104 0.0873 0.0742;...
                 1.19 0.619 0.394 0.276 0.206 0.16 0.127 0.104 0.0873 0.0742;...
                 1.19 0.619 0.394 0.266 0.177 0.13 0.101 0.0822 0.0692 0.0593];      % given in Toulorge thesis for optimized RK
Kh_DG_values = [1.45 3.42 5.56 7.78 10.1 12.4 14.7 17.1 19.4 21.8;...
                0.65 2.01 3.64 5.42 7.29 9.23 11.2 13.2 15.2 17.3;...
                0.299 1.23 2.51 3.96 5.54 7.2 8.95 10.7 12.6 14.4;...
                0.138 0.769 1.76 2.93 4.23 5.63 7.14 8.72 10.4 12.1;...
                0.0641 0.483 1.24 2.17 3.31 4.58 5.94 7.37 8.84 10.4];

% Time integration for PSTD with RKo6s
RKo6s_gamma = zeros(6,1);
RKo6s_gamma(1) = 1;
RKo6s_gamma(2) = 1/2;
RKo6s_gamma(3) = 0.165919771368;
RKo6s_gamma(4) = 0.040919732041;
RKo6s_gamma(5) = 0.007555704391;
RKo6s_gamma(6) = 0.000891421261;

% $$$ RKgamma = 1./factorial([1:Nt]);
RKo6s_coef = ones(size(RKo6s_gamma));
for ii=1:length(RKo6s_gamma)
    RKo6s_coef(end+1-ii) = RKo6s_gamma(ii)/prod(RKo6s_coef);
end

RKcoef_PS = RKo6s_coef;
RKnstage_PS = length(RKcoef_PS);
RKgamma_PS = RKo6s_gamma;

% $$$ Nt = 4;
% $$$ RKgamma = 1./factorial([1:Nt]);
% $$$ RKcoef = ones(size(RKgamma));
% $$$ for ii=1:length(RKgamma)
% $$$     RKcoef(end+1-ii) = RKgamma(ii)/prod(RKcoef);
% $$$ end