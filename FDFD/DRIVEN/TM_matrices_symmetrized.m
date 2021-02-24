function [A,b, Pr, Pl] = TM_matrices_symmetrized(wvlen, xrange, yrange, eps_r, Jz, Npml)


    %% Set up the domain parameters.
    %normal SI parameters
    eps_0 = 8.85*10^-12;
    mu_0 = 4*pi*10^-7; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
    N = size(eps_r);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% Set up the permittivity and permeability in the domain.
    % bwdmean does nearest neighbor averaging (smoothes out stuff)

    eps_z = bwdmean_w(eps0 *eps_r, 'z');
    %these are fully dense matrices...

    %currently, eps_x and eps_y are ultra-dense, which isn't right...

    %% Set up number of cells
    %the wavelength is much larger than the dimensions of the system...
    xmin = xrange(1); xmax = xrange(2);
    ymin = yrange(1); ymax = yrange(2);
    Nx = N(1); dx = (xmax-xmin)/Nx;
    Ny = N(2); dy = (ymax-ymin)/Ny;
    % Nz = 1; dz = 1; 2D solving only
    M = prod([Nx, Ny]); %total number of cells

    %% Set up the Split coordinate PML
    %sx = create_sfactor('f',Nx);
    %sy = creates_factor('f',Ny);
    Nx_pml = Npml(1); Ny_pml = Npml(2);
    Nwx = Nx; Nwy = Ny;
    sxf = create_sfactor_mine(xrange,'f',omega,eps_0,mu_0,Nwx,Nx_pml);
    syf = create_sfactor_mine(yrange,'f', omega,eps_0,mu_0,Nwy,Ny_pml);
    sxb = create_sfactor_mine(xrange, 'b', omega,eps_0,mu_0, Nwx, Nx_pml);
    syb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Nwy,Ny_pml);

    % now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
    [Sxf, Syf] = ndgrid(sxf, syf);
    [Sxb, Syb] = ndgrid(sxb, syb);

    %Sxf(:) converts from n x n t0 n^2 x 1
    Sxf=spdiags(Sxf(:),0,M,M);
    Sxb=spdiags(Sxb(:),0,M,M);
    Syf=spdiags(Syf(:),0,M,M);
    Syb=spdiags(Syb(:),0,M,M);


    %% Create the dielectric and permeability arrays (ex, ey, muz)
    %create a diagonal block matrix of ep and mu...
    epzList = reshape(eps_z,M,1);
    Tepz = spdiags(epzList,0,M,M); % creates an MxM matrix, which is the correct size,
    %the M entries in epsList is put on the diagonals
    Tmz = mu0*speye(M); %in most cases, permeability is that of free-space
    Tmy = Tmz; Tmx = Tmz;

    %% Create Magnetic vector Mz (source profile determined by Mz input)
    % dimension = M*1
    Jz = reshape(Jz,M,1);
    Jz = sparse(Jz);

    %% create the derivative oeprators w/ PML

    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws('y', 'f', dL, N);
    Dyb = createDws('y', 'b', dL, N); 
    Dxb = createDws('x', 'b', dL, N); 
    Dxf_pml = Sxf^-1*Dxf; 
    Dyf_pml = Syf^-1*Dyf;
    Dyb_pml = Syb^-1*Dyb; 
    Dxb_pml = Sxb^-1*Dxb; 


    %% Construct the matrix A, everything is in 2D
    A = Dxb_pml*(Tmy^-1)*Dxf_pml + Dyb_pml*(Tmx^-1)*Dyf_pml + omega^2*Tepz;
    % note a warning about ill-conditioned matrices will pop up here, but
    % for our purposes, it is okay.

    %% construct the matrix b, everything is in 2D
    b = 1i*omega*Jz;

    %%
    % hx = -1/(1i*omega)*(Tmx^-1*Dyf_pml)*ez;
    % hy = (Tmy^-1*Dxf_pml)*ez*(1/(1i*omega));
    hx_op = -1/(1i*omega)*(Tmx^-1*Dyf_pml);
    hy_op = (Tmy^-1*Dxf_pml)*(1/(1i*omega));

 
end