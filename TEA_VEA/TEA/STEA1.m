%% Set up the domain parameters.
close all
clear all
%% VECTOR EPSILON ALGORITHM (VEA)
%% implement a termination to prevent instability
% in evaluating large numbers of transforms.

%% ===================================================
% Here's my idea for a termination condition
% look at epsilon in epsilon_wynn, if the norm of the 
% epsilon wynn term is less than 1/max(field), then stop

%% ===================================================

wvlen = 1;
L0 = 1.55e-6;  % length unit: microns
eps0 = 8.854e-12*L0;
mu0 = 4*pi*1e-7*L0;
c0 = 1/sqrt(mu0*eps0);
Z0 = sqrt(mu0/eps0);
k0 = 2*pi/wvlen
num_cells = 2;
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
% imagesc(real(Ez_for))

II1b = 60;
II2b = 100;
                    
TTopt = ones(N);
%% ======================================================================%%
iiter = 1;
% Ez_for if the forward simulation field

% what's the different between Ez_proj and dEz_for
Ez_proj{iiter}=Ez_for;
dEz_for{1}=full(Ez_for); % dEz_for... the derivative of dEz...usualy computed by adjoint.

%II1b, II2b is the point where we are optimizing
% so are we calculating the series only for one point?
% yes, because we are only optimizing one point, so we only 
% need to calculate the field at that point for a given deltaEpsilon

S(1)=dEz_for{1}(II1b,II2b);

%% SWITCH BETWEEN N_order = 20
%% N_order = 7 should produce convergent...
N_order=21; %13;

%% deepstr
deeps_str = zeros(N);
%deeps_str(40:60, 40:60) = 3;
% deeps_str(50:90, 50:90) = 2;
%deeps_str(60:80, 60:80) = 4;
deeps_str(40:100, 40:100) = 1;
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
    dEz_for{iorder} = full(new_field);
    %[dEz_for{iorder}] = solveTM(wvlen, xrange, yrange, eps_str_old, Mz, Npml);

end;
toc

figure();
imagesc(real(Mz_LS_cell{2})); colorbar;

% I suspect this is the search grid, but not sure
Shanks_cell{1} = dEz_for{1};

%% correctly, you don't have to re-evaluate the green's function beyond this point

alphav = [0.5, 0.8]
for ialpha = 1:length(alphav)

%% for lippman-schwinger??
for iorder=2:N_order
    %% lipmman schwinger prediction of the field
    Ez_proj{iiter}=Ez_proj{iiter}+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder};

    %% initialization of the shanks transform of the field
    S(iorder)=S(iorder-1)+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder}(II1b,II2b);
    Shanks_cell{iorder} = Shanks_cell{iorder-1} + (alphav(ialpha)).^(iorder-1)*dEz_for{iorder};
end


%% for shanks??
e{1}=S; %stored as a cell, where cell consists of ip and then each cell contains a value
        % it would seem
eShanks{1} = Shanks_cell;
eShanks{2} = {};

%% ====================== DEFINE OUR y vect for bilinear ================%%
%remarkable degree of insensitivity to y
y = full(Ez_for);
%y  = abs(10*randn(N));
%y(zone==0) = 0;

%% calculate the SEA of the bilinear form <y, Sn>
%% this is the most expensive part of the algorithm...but why?
tic;
%% why is this so expensive?
yS = zeros(N_order,1);
for i = 1:N_order
    %% norming is actually a relatively expensive operation
    % it is two operations, squaring and modulus ds
    yS(i) = norm(y.*eShanks{1}{i});
end
toc;
e{1}=yS;
for iorder = 1:N_order-1
    e{2}(iorder)=1/(yS(iorder+1)-yS(iorder));
end
for ip=3:N_order
    for iorder=1:N_order-ip+1
        % this looks much like the shanks transformation
        denominator = 1/(e{ip-1}(iorder+1)-e{ip-1}(iorder));
        e{ip}(iorder)=e{ip-2}(iorder+1)+1/(e{ip-1}(iorder+1)-e{ip-1}(iorder));        
    end
end

%% STEA1 requires only even 2k orders?
% it is at least twice as fast as the original TEA1, but it is still too
% slow
tic
stea{1} = Shanks_cell;
stea{2} = {};
for ip=3:2:N_order
    for iorder=1:N_order-ip+1
        scalar_term = 1/((e{ip-2}(iorder+1)-e{ip-2}(iorder))*(e{ip-1}(iorder+1)-e{ip-1}(iorder)));
        stea{ip}{iorder} = stea{ip-2}{iorder+1}+(scalar_term)*(stea{ip-2}{iorder+1}-stea{ip-2}{iorder});
    end
end
toc
figure(); 
imagesc(real(stea{N_order}{1}))

%% ORIGINAL TEA1
tic; %% this is the cheaper step
for iorder=1:N_order-1
    delta1n = eShanks{1}{iorder+1}-eShanks{1}{iorder};
    %samelson_inverse = conj(z)/norm(full(z))^2;
    denominator = conj(delta1n).*y;
    teaInverse = y./denominator;    
    eShanks{2}{iorder} = teaInverse;
end
toc
tic
%% now starting from n=2
denominator_norms = [];
for ip=3:N_order
    eShanks{ip} = {};
    for iorder=1:N_order-ip+1
        
        %% for tea, we have to seperate out even and odd in ip or iorder
        if(mod(ip-1,2) == 0)%2k+2
            delta1n = eShanks{ip-1}{iorder+1}-eShanks{ip-1}{iorder};
            delta0n = eShanks{ip-2}{iorder+1}-eShanks{ip-2}{iorder};
            bilinear = conj(delta1n).*delta0n; %THIS ORDER IS IMPORTANT
            teaInverse = delta0n./bilinear;
        else %2k+1, odd cases
            delta1n = eShanks{ip-1}{iorder+1}-eShanks{ip-1}{iorder};
            %samelson_inverse = conj(z)/norm(full(z))^2;
            denominator = conj(delta1n).*y;            
            teaInverse = y./denominator;
        end
        eShanks{ip}{iorder} = eShanks{ip-2}{iorder+1} + teaInverse;

    end

end
toc;

%% validate prediction
new_eps = alphav(ialpha)*deeps_str+eps_str_old;
Ez_validation = solveTM(L0, wvlen, xrange, yrange, new_eps, Mz0, Npml);

figure();
subplot(1,3,1)
imagesc(real(Ez_validation)); colorbar

%visreal(Ez_validation, xrange, yrange); colorbar
title(strcat('alpha=',num2str(alphav(ialpha))));
valcaxis = caxis;
subplot(1,3,2)
imagesc(real(eShanks{N_order}{1})); colorbar;

%visreal(real(eShanks{N_order}{1}), xrange, yrange); colorbar;
caxis(valcaxis)
%figure();
subplot(1,3,3)
imagesc(real(stea{N_order}{1})); colorbar;
title('stea1')
%visabs(abs(eShanks{N_order}{1}-Ez_validation)/(max(max(abs(Ez_validation)))), xrange, yrange); colorbar;


%% check the actual series
%% Later shanks cell terms are LARGER in magnitude (for faulty cases)
% is that correct? NO...good cases tend to CONTROL the growth 
% in later field terms, but why?
% for i = 1:N_order
%    figure()
%    imagesc(real(Shanks_cell{i}))
%    colorbar();
%    drawnow();
% end

%% check the shanks transform that feeds into N_order
% for i = 3:2:N_order
%     figure();
%     imagesc(real(eShanks{i-2}{end}));
%     title(i) 
% end

%% numerical errors are already appearing
figure();
subplot(1,2,1)
imagesc(real(eShanks{3}{3})); colorbar;
subplot(1,2,2)
imagesc(real(eShanks{5}{1})); colorbar;

%% terms 1 and 5 have artifacts...
figure();
subplot(1,2,1)
imagesc(real(eShanks{3}{1}))
subplot(1,2,2)
imagesc(real(eShanks{3}{3}))

%
%% eShanks{3}{5} = eShanks{1}{5+1} + 1./(eShanks{2}{5+1} - eShanks{2}{5});
% first term is usually okay
figure();
imagesc(real(eShanks{1}{6}));

%% numerical stabilities are sourced from the 
%% even terms in this array
figure();
subplot(1,2,1)
imagesc(real(1./(eShanks{2}{5+1} - eShanks{2}{5})));
colorbar();
subplot(1,2,2)
imagesc(real((eShanks{2}{5+1} - eShanks{2}{5})));
colorbar();

figure();
subplot(1,2,1)
imagesc(log10(abs(eShanks{2}{5+1} ))); %% this one has some error
colorbar();
subplot(1,2,2)
%% the pml causes numerical problems...
imagesc(log10(abs((eShanks{2}{5}))));
colorbar();

%% eShanks{2}{iorder} = 1./(Shanks_cell{iorder+1} - Shanks_cell{iorder});
figure();
subplot(1,2,1)
imagesc(real(Shanks_cell{5+1}))
subplot(1,2,2)
imagesc(real(Shanks_cell{5}))

%% probably want to stabilize numerical
figure();
subplot(1,2,1)
imagesc(real((Shanks_cell{5+1} - Shanks_cell{5})))
subplot(1,2,2)
imagesc(log10(1./abs((Shanks_cell{5+1} - Shanks_cell{5}))));

colorbar();

end