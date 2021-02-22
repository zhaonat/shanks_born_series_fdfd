%% Set up the domain parameters.

%% numerator/denominator implementation of the Shanks transformation

close all
clear all

wvlen = 1;
L0 = 1.55e-6;  % length unit: microns
eps0 = 8.854e-12*L0;
mu0 = 4*pi*1e-7*L0;
c0 = 1/sqrt(mu0*eps0);
Z0 = sqrt(mu0/eps0);
k0 = 2*pi/wvlen
xrange = [-1.25 1.25];  % x boundaries in L0
yrange = [-1.25 1.25];  % y boundaries in L0

L = [diff(xrange), diff(yrange)];
N = [125 125];  % [Nx Ny]
N0 = N; 
Npml = 1*[10 10];  % [Nx_pml Ny_pml]

zone = zeros(N);
zone(Npml(1):end-Npml(1),Npml(2):end-Npml(2)) = 1;
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
[X,Y]=meshgrid(linspace(xrange(1),xrange(2),N(1)),linspace(yrange(1),yrange(2),N(2)));


eps_r = ones(N);
eps_str_old = eps_r;

Mz = zeros(N);
Mz(20,70) = 1;

[A,b, omega] = solveTM_matrices(wvlen, xrange, yrange, eps_str_old, Mz, Npml);

tic
[L,U,P,Q] = lu(A);
toc

tic
Ez_for = reshape(P.'*(U\(L\(Q.'*b))),N);
toc
Mz0 = Mz;


II1b = 60;
II2b = 100;
                    
TTopt = ones(N);
%% ======================================================================%%

% what's the different between Ez_proj and dEz_for
Ez_proj=Ez_for;
dEz_for{1}=Ez_for; % dEz_for... the derivative of dEz...usualy computed by adjoint.

%II1b, II2b is the point where we are optimizing
% so are we calculating the series only for one point?
% yes, because we are only optimizing one point, so we only 
% need to calculate the field at that point for a given deltaEpsilon

S(1)=dEz_for{1}(II1b,II2b);

%% SWITCH BETWEEN N_order = 20
%% N_order = 7 should produce convergent...
N_order=7;

%% deepstr
deeps_str = zeros(N);
%deeps_str(40:60, 40:60) = 3;
deeps_str(50:90, 50:90) = 1; %randn(41,41);
%deeps_str(60:80, 60:80) = 4;
%deeps_str(40:100, 40:100) = 1;
%deeps_str(10:end-10, 10:end-10)= 1;
%% it appears we are doing several forward simulations here to get the sources%deeps_str(40:60, 40:60) = 3;

% at differenert iorders, but why?
Mz_LS_cell{1} = Mz0;
tic
for iorder=2:N_order
   
    % even if I set TTopt to 1, the code below extracts the region where 
    % deeps_str > 0, i.e. it extracts only the region of the perturbation
    % to construct a source
    
    % only setting source in the Topt region?
    Mz(TTopt==1)=1./(-1j*k0*Z0)*k0^2*deeps_str(TTopt==1).*dEz_for{iorder-1}(TTopt==1);
    Mz_LS_cell{iorder} = Mz;
    bprime = sparse(reshape(Mz, prod(N), 1));
    bprime = (1j*omega)*bprime;
    
    new_field = reshape(P.'*(U\(L\(Q.'* bprime))),N);
    dEz_for{iorder} = new_field;
    %[dEz_for{iorder}] = solveTM(wvlen, xrange, yrange, eps_str_old, Mz, Npml);

end;
toc

figure();
imagesc(real(Mz_LS_cell{2})); colorbar;

% I suspect this is the search grid, but not sure
Shanks_cell{1} = dEz_for{1};

alphav = [0.8]
ialpha = 1

%% for lippman-schwinger??

for iorder=2:N_order
    %% lipmman schwinger prediction of the field
    Ez_proj=Ez_proj+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder};

    %% initialization of the shanks transform of the field
    S(iorder)=S(iorder-1)+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder}(II1b,II2b);
    Shanks_cell{iorder} = Shanks_cell{iorder-1} + (alphav(ialpha)).^(iorder-1)*dEz_for{iorder};

end


%% SHANKS transform
shanks_order = 2;
eShanks=cell(1);
eShanks{1} = Shanks_cell;

for sorder = 1:shanks_order
    n = length(eShanks{sorder});
    for in = 2:n-1
        numerator = eShanks{sorder}{in+1}.*eShanks{sorder}{in-1} -(eShanks{sorder}{in}).^2;
        denominator = eShanks{sorder}{in+1} + eShanks{sorder}{in-1}-2*eShanks{sorder}{in};
        eShanks{sorder+1}{in-1} = numerator./denominator;
    end

end


%% validate prediction
new_eps = deeps_str+eps_str_old;
Ez_validation = solveTM(wvlen, xrange, yrange, new_eps, Mz0, Npml);


figure();
imagesc(real(Ez_validation));
colorbar();
valcaxis = caxis

figure();
imagesc(real(eShanks{3}{2}));
colorbar;
caxis(valcaxis);
