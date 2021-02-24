function [A, b, omega, Sxf, Syf, Sxb, Syb, Dxf, Dyf, T_eps_x, T_eps_y] = ...
    solveTE_Jdipole_matrices(L0, wvlen, xrange, yrange, eps_r, Mz, Npml, Jx, Jy)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Mz: Nx-by-Ny array of magnetic current source density
% Jx, Jy: Nx-by-Ny array of ELECTRIC current source
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters
% Ez, Hx, Hy: Nx-by-Ny arrays of H- and E-field components
% dL: [dx dy] in L0
% A: system matrix of A x = b
% omega: angular frequency for given wvlen

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./N;  % [dx dy]

M = prod(N); 

omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec

%% Deal with the s_factor

%% Set up the permittivity and permeability in the domain.
eps_x = bwdmean_w(eps0*eps_r, 'x'); 
eps_y = bwdmean_w(eps0*eps_r, 'y'); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
T_eps_x = spdiags(eps_x(:), 0, M, M); 
T_eps_y = spdiags(eps_y(:), 0, M, M); 


%% Construct derivate matrices
Dyb = Syb\createDws('y', 'b', dL, N); 
Dxb = Sxb\createDws('x', 'b', dL, N); 
Dxf = Sxf\createDws('x', 'f', dL, N); 
Dyf = Syf\createDws('y', 'f', dL, N); 

%% Reshape Mz into a vector
mz = reshape(Mz, M, 1); 
jx = reshape(Jx, M, 1);
jy = reshape(Jy, M, 1);

%% contribution of j
jx_con = -Dyf*(T_eps_x^-1)*jx;
jy_con = Dxf*(T_eps_y^-1)*jy;

%% Construct A matrix and b vector
% A = Sx_f *Dxf* T_eps_y^-1 *Sx_b *Dxb + Sy_f *Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z; 
% A = Sx_f*Dxf* T_eps_y^-1 *Sx_b*Dxb + Sy_f*Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z; 

I = speye(M); 
A_mode = Dxf*(T_eps_x^-1)*Dxb + Dyf*(T_eps_y^-1)*Dyb;
A = A_mode + omega^2*mu0*I; 
% % A = Dxf* T_eps_y^-1 *Dxb + Dyf* T_eps_x^-1* Dyb + omega^2*T_mu_z; 
b = 1i * omega * mz + jx_con + jy_con; 


end
