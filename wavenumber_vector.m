function kxvec = wavenumber_vector(Nx, dx)
% WAVENUMBER_MAT: build discrete wavenumber vector for PSTD computations
% Nx: number of points along the considered axis
% dx: sampling interval

kxvec = (0:Nx-1)/(Nx*dx);
kxvec(kxvec>1/2/dx) = kxvec(kxvec>1/2/dx)-1/dx;
kxvec = 2*pi*kxvec;
