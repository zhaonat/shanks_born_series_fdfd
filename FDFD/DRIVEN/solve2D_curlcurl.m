function [Ex, Ey, Ez, Hx, Hy, Hz, A, omega,b, Sxf, Dxf, Dyf, sxf, syf,trun] = ...
    solve2D_curlcurl(wvlen, xrange, yrange, eps_r, JCurrentVector, Npml)
    % SOLUTION  using full 3D curl curl operator... just for demonstration
    % don't actually use this for 2D simulations.
    
    %% Input Parameters
    % wvlen: wavelength in L0
    % xrange: [xmin xmax], range of domain in x-direction including PML
    % yrange: [ymin ymax], range of domain in y-direction including PML
    % eps_r: Nx-by-Ny array of relative permittivity
    % Mz: Nx-by-Ny array of magnetic current source density
    % Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

    %% Output Parameters
    % Hz, Ex, Ey: Nx-by-Ny arrays of H- and E-field components
    % dL: [dx dy] in L0
    % A: system matrix of A x = b
    % omega: angular frequency for given wvlen

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

    eps_x = bwdmean_w(eps0 * eps_r, 'y');  % average eps for eps_x
    eps_y = bwdmean_w(eps0 * eps_r, 'x');  % average eps for eps_y
    eps_z = bwdmean_w(eps0 * eps_r, 'z');
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
    epxList = reshape(eps_x,M,1);
    epyList = reshape(eps_y,M,1);
    epzList = reshape(eps_z,M,1);
    Tepx = spdiags(epxList,0,M,M); % creates an MxM matrix, which is the correct size,
    %the M entries in epsList is put on the diagonals
    Tepy = spdiags(epyList,0,M,M);
    Tepz = spdiags(epzList, 0,M,M);
    Tmz = mu0*speye(M); %in most cases, permeability is that of free-space
    TepsSuper = blkdiag(Tepx,Tepy,Tepz);
    TmuSuper = blkdiag(Tmz, Tmz, Tmz);
    %% Create Magnetic vector Mz (source profile determined by Mz input)
    % dimension = M*1
    Mx = JCurrentVector(1:Nx,:,:); My = JCurrentVector(Nx+1:2*Nx,:,:); Mz = JCurrentVector(2*Nx+1: 3*Nx,:,:);
    Mz = reshape(Mz, M, 1);
    My = reshape(My, M, 1);
    Mx = reshape(Mx, M, 1);
    Mz = sparse(Mz);
    J = [Mx; My; Mz];
   
    %% create the derivative oeprators w/ PML

    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws_dense('x', 'f', dL, N); 
    Dyf = createDws_dense('y', 'f', dL, N);
    Dyb = createDws_dense('y', 'b', dL, N); 
    Dxb = createDws_dense('x', 'b', dL, N); 
    Dxf= Sxf^-1*Dxf; 
    Dyf = Syf^-1*Dyf;
    Dyb = Syb^-1*Dyb; 
    Dxb= Sxb^-1*Dxb; 
    
    %% Construct Wonsoek's Accelerator
    GradDiv = [Dxf_pml*Dxf_pml Dxf_pml*Dyf_pml zeros(M); Dyf_pml*Dxf_pml Dyf_pml*Dyf_pml zeros(M); ...
        zeros(M) zeros(M) zeros(M)];
    
    GradDiv = sparse(GradDiv);
    
    Ce = [zeros(M) zeros(M) Dyf; zeros(M) zeros(M) -Dxf; -Dyf Dxf zeros(M)];
    Ch = [zeros(M) zeros(M) Dyb; zeros(M) zeros(M) -Dxb; -Dyb Dxb zeros(M)];
    Ce = sparse(Ce);
    Ch = sparse(Ch);
    
    %% PROBLEM CONSTRUCTION
    b = -1i*omega*J;
    JCorrection = s*(1i/omega) * GradDiv*J;
    b = b+JCorrection;
    A = Ch*TmuSuper^-1*Ce - s*GradDiv - omega^2*TepsSuper;
    solution = A\b;
    
    %% SOLUTION EXTRAction
    solLength = length(solution);
    Ex = solution(1:solLength/3);
    Ey = solution(solLength/3+1:solLength*(2/3));
    Ez = solution(solLength*(2/3)+1: solLength);
    h = (1i/omega)*TmuSuper^-1*(Ce*solution);
    
    Hx = h(1:solLength/3);
    Hy = h(solLength/3+1:solLength*(2/3));
    Hz = h(solLength*(2/3)+1: solLength);

 
end