close all
clear all

%% gradient scaling or bound constraints are critical
% gradient clipping effectively destroys any benefit the
% shanks transformation might give you

%% Set up the domain parameters.
L0 = 1.55e-6;  % length unit: microns
eps0 = 8.854e-12*L0;
mu0 = 4*pi*1e-7*L0;
Z0 = sqrt(mu0/eps0);
c0 = 1/sqrt(mu0*eps0);                     
xrange = [-1.25 1.25];  % x boundaries in L0
yrange = [-1.25 1.25];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [125 125];  % [Nx Ny]
N0 = N; 
Npml = 1*[10 10];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
[X,Y]=meshgrid(linspace(xrange(1),xrange(2),N(1)),linspace(yrange(1),yrange(2),N(2)));


Nx = N(1); Ny = N(2);
M = prod(N);
%% set up initial epsilon
eps_r= ones(N).*1.5^2;

%% source placement
wvlen = 1;
k0 = 2*pi/wvlen
omega = 2*pi*c0/wvlen;
source_coord = [0,0];
[njx, njy] = coord_to_grid(source_coord, N, xrange, yrange);
Jz = zeros(N);
Jz(njx, njy) = 1;

%% define coordinate to optimize
opt_coord = [0, 1];
[nox, noy] = coord_to_grid(opt_coord, N, xrange, yrange);
eta = zeros(N); 
eta(nox, noy) = 1;
%eta = eta(:);

%% II1b, II2b is Salim's notation for the integer coordinates.
II1b = nox;
II2b = noy;

%% ============== DEFINE Schur complement====================%%
%% define a subdomain rectangle
corner1 = [-0.5, -0.5];
corner2 = [0.5, 0.5];
measure_point = [0,0];

[cx1, cy1] = coord_to_grid(corner1, N, xrange, yrange);
[cx2, cy2] = coord_to_grid(corner2, N, xrange, yrange);
coord1 = [cx1, cy1];
coord2 = [cx2, cy2];
% the optimization subdomain is different, defined by corner1, corner2
opt_size = coord2-coord1-1;
xrange_sub = [-0.5, 0.5];
yrange_sub = [-0.5, 0.5];
[nxs, nys] = coord_to_grid(measure_point, opt_size, ...
    xrange_sub, yrange_sub);


%% ENTER Optimization loop
epochs = 1000;%250;%1000;%NB Salim number of iterations

alpha = 1;%NB permittivity step factor

%% define region inside and outside PML
zone=zeros(N); 
zone(((X>=-1).*(X<=1).*(Y>=-1).*(Y<=1))==1)=1;
zone=reshape(zone,N);


%% ========================= MOMENTUM TERM ===============================
vdw = zeros(N); %% this stores our momentum
beta = 0.65;

%%  ========================= RMSProp TERM ===============================
past_derivative_of_objective = ones(N);
increment_factor = 1.1;
decrement_factor = 0.9;
max_alpha = 0.4;
alpha_grid = max_alpha*ones(N);
alpha_grid(zone == 0) = 0;

%%  =========================== ADAM TERMS =============================
epochs = 100;
shanks_alpha_history = [];
prev_field = 0;
for t =1:epochs
    tic
    t0=cputime;
    [A,b, omega] = ...
    solveTM_matrices(wvlen, xrange, yrange, eps_r, Jz, Npml);
    [L,U,P,Q] = lu(A);
    Ez =  A\b;
    Ez=reshape(Ez,N);
    prev_field = Ez;

    if t==1
        ref=log10((abs(Ez(eta==1)))^2);FOM(t)=ref;
    else,
        FOM(t) = log10((abs(Ez(eta==1)))^2);
    end;
%     FOM(t) = conj(Ez(eta==1))*(Ez(eta==1));
    grad_objective = reshape(eta.*conj(Ez),M,1);%-2*(eta.'*u0)*eta;%NB Salim the adjoint formulation was wrong. I corrected it.

    Ez_adj = A\grad_objective;
    Ez_adj=reshape(Ez_adj,N);

    derivative_of_objective = full(real(Ez.*Ez_adj));%FOM(t-1)*100*max(max(abs(full(omega^2*(eps0/L0)*real(Ez.*Ez_adj))))));%full(-omega^2*(eps0/L0)*real(Ez.*Ezu));%NB Salim the adjoint formulation was wrong.
    d0 = derivative_of_objective;
    %% what does this scaling of the derivative of the objective
    derivative_of_objective =derivative_of_objective./(max(abs(derivative_of_objective)));
    
    deeps_str = -derivative_of_objective;
    
    if(t==1)
        Time(t)=cputime-t0;
    else
        Time(t)=Time(t-1)+cputime-t0;
    end
    toc
    
    %% add in line search 
    alpha_scan = logspace(-2, -0.01, 20);
    N_order = 5;
    best_new_fom = FOM(1);
    best_alpha = alpha_scan(1);
    best_field = Ez;
    
    masked_deeps_str = deeps_str;
    masked_deeps_str(zone == 0) = 0;
    alphav = alpha_scan;
    alphai = 1;
    
    [dEz_for] = lippman_schwinger(Ez, L,U,P,Q, Jz, masked_deeps_str, N_order, k0, Z0, omega);
    figure()
    for alpha = alpha_scan
        eShanks =...
          shanks_transformation(dEz_for, alpha, N_order);
      new_field = eShanks{N_order}{1};
      fom = log10((abs(new_field(eta==1)))^2);
      prediction(alphai) = fom;
      if(fom > best_new_fom)
         best_new_fom = fom;
         best_alpha = alpha;
         best_field = new_field;
      end
%       if( mod(alphai,10) == 0)
% 
%          subplot(1,4,1)
%          title('original')
%          imagesc(abs(Ez))
%          cvalaxis = caxis;
%          colorbar;
%          subplot(1,4,2)
%          imagesc(abs(new_field))
%          title('new')
%          caxis(cvalaxis); colorbar;
%          subplot(1,4,3);
%          plot(alpha_scan(1:length(prediction)),prediction)
%          subplot(1,4,4);
%          imagesc(alpha*derivative_of_objective);
%          drawnow()
%       end
      alphai = alphai+1;
    end
    shanks_alpha_history = [shanks_alpha_history, best_alpha];
    prev_field = best_field;
    %% ==========================================================%%
    %% simple gradient;
    dw = best_alpha*derivative_of_objective;
    dw(zone==0) = 0;
    %% momentum...term...
    vdw(zone>0) = beta*vdw(zone>0) + alpha*derivative_of_objective(zone>0);
    
    %% nesterov momentum (computes gradient at "updated" new position...
    %% not sure we can do nesterov here...
    
    %% rmsprop (increment or decrement alpha based on the previous history of alphas)
    signs = -1+2*((derivative_of_objective.*past_derivative_of_objective)> 0);
    step_increment = zeros(N);
    step_increment(signs>0) = increment_factor;
    step_increment(signs<0) = decrement_factor;
    alpha_grid(zone>0) = step_increment(zone>0).*alpha_grid(zone>0);
    alpha_grid(alpha_grid > max_alpha) = max_alpha;
    alpha_grid(alpha_grid < -max_alpha) = -max_alpha;
    
    %% smoothing test of alpha grid (does not work)
    % alpha_grid(zone>0) = imgaussfilt(alpha_grid(zone>0),2);
    
    update_rms = zeros(N);
    update_rms(zone>0) = alpha_grid(zone>0).*derivative_of_objective(zone>0);
    past_derivative_of_objective = derivative_of_objective;
    
    %% adam THE MOST ADVANCED POSSIBLE ALGORITHM
    %%=================== PICK UPDATE ==================$$
    update = update_rms;
    %% only update the zone of IMAX...with this simulation.
%     new_eps = alpha*derivative_of_objective(zone==Imax);
%     eps_r(zone == Imax)=eps_r(zone==Imax)- new_eps;
%     old_eps_r = eps_r;
%     eps_r = eps_r - update;
%     eps_r(eps_r<1.5^2)=1.5^2;
%     eps_r(eps_r>2.5^2)=2.5^2;
%     if(mod(t,10) == 0)
%        alpha = alpha*0.95; 
%     end
    


    if(t>1 && FOM(t) < FOM(t-1))
       alpha = alpha*0.5; 
    end
    Ez_for = Ez;
    eps_str = eps_r;
    TTopt = zone;
    eps_r = eps_r - dw;
    
%     %% =============== SPECIALIZED UPDATE =================%%
    %% this is a very different way to compute FOM...
    dFOM=real((2.5^2-1.5^2).*Ez_for.*conj(Ez_for(II1b,II2b)).*Ez_adj);
    dFOM(dFOM>0)=dFOM(dFOM>0).*(2.5^2-eps_str(dFOM>0))./((2.5^2-1.5^2));
    dFOM(dFOM<0)=dFOM(dFOM<0).*(eps_str(dFOM<0)-1.5^2)./((2.5^2-1.5^2));

%     %% this is to maintain constraint... but I do not see it used.
% %     term1 = (2.5^2-eps_str((TTopt.*(dFOM>0).*(eps_str>=(1.5^2))...
% %         .*(eps_str<=(2.5^2)))==1));
% %     term2 = dFOM((TTopt.*(dFOM>0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1);
% %     
% %     term3 = eps_str((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)-1.5^2)./
% %     term4 = (-dFOM((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)))
% % %     
% %     facmin=min(min((2.5^2-eps_str((TTopt.*(dFOM>0).*(eps_str>(1.5^2))...
% %         .*(eps_str<(2.5^2)))==1))./dFOM((TTopt.*(dFOM>0)...
% %         .*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)), ...
% %         min((eps_str((TTopt.*(dFOM<0).*(eps_str>(1.5^2))...
% %         .*(eps_str<(2.5^2)))==1)-1.5^2)./(-dFOM((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1))));
%     
%     facmin=min(min((2.5^2-eps_str((TTopt.*(dFOM>0).*(eps_str>=(1.5^2))...
%         .*(eps_str<=(2.5^2)))==1))./dFOM((TTopt.*(dFOM>0)...
%         .*(eps_str>=(1.5^2)).*(eps_str<=(2.5^2)))==1)), ...
%         min((eps_str((TTopt.*(dFOM<0).*(eps_str>=(1.5^2))...
%         .*(eps_str<=(2.5^2)))==1)-1.5^2)./(-dFOM((TTopt.*(dFOM<0)...
%         .*(eps_str>=(1.5^2)).*(eps_str<=(2.5^2)))==1))));
%     
% 
%     %% interpolation scheme for the choice of best alpha
%     UUUint=1e-6:1e-6:1;
%     VVVint=interp1(alphav,full(prediction),UUUint,'spline');
%     %% this is the best alphamax based on the interpolation
%     WWW=UUUint(find(VVVint==max(VVVint)));
%     alphamax=WWW;
% 
%     %the best change in permittivity
%     eps_str(TTopt==1)=eps_str(TTopt==1)+dFOM(TTopt==1).*(alphamax).*facmin;
%     %eps_r = eps_str;

    %% ==============================================
    
%     if(mod(t, 5) == 0)
%         figure(2);
%         pcolor(X,Y,eps_r.')
%         title('Epsilon')
%         shading('interp')
%         colormap(flipud(colormap('gray')));
%         colorbar
%         drawnow
% 
%         figure(3);
%         pcolor(X,Y,abs(Ez.'))
%         title('Ez')
%         shading('interp')
%         colorbar
%         drawnow
% 
%         figure(4);
%         plot(FOM)
%         xlabel('Iteration number')
%         ylabel('log_{10}(FOM)')
%         drawnow
% 
%         figure(5);
%         subplot(121)
%         plot(Time,FOM-ref)
%         xlabel('Cumputation time (s)')
%         ylabel('log_{10}(FOM)')
%         subplot(122)
%         plot(Time,10.^(FOM-ref))
%         xlabel('Cumputation time (s)')
%         ylabel('(FOM)')
%         drawnow
%     end
    if(t>1)
       figure();
       subplot(1,4,1);
       imagesc(real(prev_field));
       title('Shanks')
       colorbar();
       subplot(1,4,2);
       imagesc(real(Ez));
       colorbar();
       title('actual')
       subplot(1,4, 3);
       imagesc(abs(Ez-prev_field)); colorbar;
     
       subplot(1,4,4)
       imagesc(eps_r); colorbar();
       drawnow();
    end
    
end

%% reconstruct field at end
%recField = reconstruct_four_block(u0, Avv, Avp, bv, Q);
%Hz = reshape(recField, N(1), N(2));

%% test the shanks field
tic;
[L,U,P,Q] = lu(A);
toc;
