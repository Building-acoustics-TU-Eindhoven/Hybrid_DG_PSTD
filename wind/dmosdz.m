% Monin-Obukhov Similarity relationships (thanks to Adrien)
% Based on Cheinet, 2012
function [duxdz] = dmosdz(z, ustar, Fb, z0, T0, zref, Tref, phi)

    global g kappa
    g = 9.81; % constant of gravitation
    kappa = 0.4; % Von Karman constant

    % ustar: frictional velocity (Meppen, [0.03; 0.3])
    % Fb:    buoyancy flux (Meppen, [-0.01; 0.04]
    % z0:    roughness length (Meppen, z0 = 0.005)
    % T0:    reference surface temperature
    % zref:  zref : reference height
    % Tref:  temperature at zref

    z = z + z0;
    %z = max(z,z0);

    %[u_av] = compute_MOS(z, ustar, Fb, z0, T0, zref, Tref);
    %[dudz_av] = compute_dMOSdz(abs(z), ustar, Fb, z0, T0, zref, Tref);
    [dudz_av, d2udz2_av] = compute_dMOSdz_sym(abs(z), ustar, Fb, z0, T0, zref, Tref);

    %ux = u_av * cos(phi);
    duxdz = dudz_av * cos(phi);
    %d2uxdz2 = d2udz2_av * cos(phi);

end


function val = phi_m(zeta,L)
    if L>0
        val = -5*zeta;
    elseif L<0
        tmp = (1-16*zeta).^(0.25);
        val = 2*log(.5*(1+tmp)) + log(.5*(1+tmp.^2)) - ...
              2*atan(tmp) + pi/2;
    else
        val = 0;
    end
end


function val = phi_h(zeta,L)
    if L>0
        val = -5*zeta;
    elseif L<0
        tmp = (1-16*zeta).^(0.5);
        val = 2*log(.5*(1+tmp));
    else
        val = 0;
    end
end



function [u_av] = compute_MOS(z, ustar, Fb, z0, T0, zref, Tref)
% z    : height
% ustar: frictional velocity
% Fb:    buoyancy flux
% z0:    roughness length
% T0: reference surface temperature
% zref : reference height
% Tref : temperature at zref

    global g kappa

    Tstar = -Fb/ustar; % frictional temperature
    L = -ustar^3*T0/kappa/g/Fb; % Monin-Obukhov length
    Gamma_d = -0.0098; % dry adiabatic lapse rate

    u_av = ustar/kappa*(log(z/z0)-phi_m(z/L,L)+phi_m(z0/L,L));

end


% $$$ function [dudz_av] = compute_dMOSdz(z, ustar, Fb, z0, T0, zref, Tref)
% $$$ % z    : height
% $$$ % ustar: frictional velocity
% $$$ % Fb:    buoyancy flux
% $$$ % z0:    roughness length
% $$$ % T0: reference surface temperature
% $$$ % zref : reference height
% $$$ % Tref : temperature at zref
% $$$
% $$$     global g kappa
% $$$
% $$$     Tstar = -Fb/ustar; % frictional temperature
% $$$     L = -ustar^3*T0/kappa/g/Fb; % Monin-Obukhov length
% $$$     Gamma_d = -0.0098; % dry adiabatic lapse rate
% $$$
% $$$     dudz_av = ustar/kappa*(1./z-dphi_mdz(z/L,L));
% $$$
% $$$ end


    function [dudz_av, d2udz2_av] = compute_dMOSdz_sym(z, ustar, Fb, z0, T0, zref, Tref)
    % z    : height
    % ustar: frictional velocity
    % Fb:    buoyancy flux
    % z0:    roughness length
    % T0: reference surface temperature
    % zref : reference height
    % Tref : temperature at zref

        global g kappa

        Tstar = -Fb/ustar; % frictional temperature
        L = -ustar^3*T0/kappa/g/Fb; % Monin-Obukhov length
        Gamma_d = -0.0098; % dry adiabatic lapse rate

        for ii = 1:size(z,1)
            for jj = 1:size(z,2)

                zii = z(ii,jj);

                % first derivative
                if L<0
                    dudz_av(ii,jj) = (ustar*(1/zii - 8/(L*((1 - (16*zii)/L)^(1/2) + 1)*(1 ...
                                                                      - ...
                                                                      (16*zii)/L)^(3/4)) ...
                                             + 4/(L*((1 - (16*zii)/L)^(1/2)/2 + 1/2)*(1 ...
                                                                      - ...
                                                                      (16*zii)/L)^(1/2)) ...
                                             + 4/(L*((1 - (16*zii)/L)^(1/4)/2 + 1/2)*(1 ...
                                                                      - ...
                                                                      (16*zii)/ ...
                                                                      L)^(3/4))))/kappa;
                elseif L>=0
                    dudz_av(ii,jj) = (ustar*(5/L + 1/zii))/kappa;
                end

                % second derivative
                if L<0
                    d2udz2_av(ii,jj) = -(ustar*(1/zii^2 + 64/(L^2*((1 - (16*zii)/L)^(1/2) + 1)^2*(1 - ...
                                                                      (16*zii)/L)^(5/4)) ...
                                                + 96/(L^2*((1 - (16*zii)/L)^(1/2) + 1)*(1 - (16*zii)/ ...
                                                                      L)^(7/4)) + ...
                                                16/(L^2*((1 - (16*zii)/L)^(1/2)/2 + 1/2)^2*((16*zii)/L ...
                                                                      - 1)) - ...
                                                32/(L^2*((1 - (16*zii)/L)^(1/2)/2 + 1/2)*(1 - ...
                                                                      (16*zii)/L)^(3/2)) ...
                                                - 8/(L^2*((1 - (16*zii)/L)^(1/4)/2 + 1/2)^2*(1 - ...
                                                                      (16*zii)/L)^(3/2)) ...
                                                - 48/(L^2*((1 - (16*zii)/L)^(1/4)/2 + 1/2)*(1 - ...
                                                                      (16*zii)/ ...
                                                                      L)^(7/ ...
                                                                      4))))/kappa;
                elseif L>=0
                    d2udz2_av(ii,jj) = -ustar/(kappa*zii^2);
                end
            end
        end

    end
