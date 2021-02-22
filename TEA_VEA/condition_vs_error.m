%% Set up the domain parameters.
close all
clear all


%% estimated condition number is ~1e5 for this problem.
% so five digits of accuracy lost per factorization?

%% numerical artifacts get worse with smaller wavelengths
wvlen_scan = linspace(0.5, 1.5, 40);
for wi = 1:length(wvlen_scan)
    close all
    wvlen = wvlen_scan(wi);
    L0 = 1.55e-6;  % length unit: microns
    eps0 = 8.854e-12*L0;
    mu0 = 4*pi*1e-7*L0;
    c0 = 1/sqrt(mu0*eps0);
    Z0 = sqrt(mu0/eps0);
    k0 = 2*pi/wvlen;
    xrange = [-1.25 1.25];  % x boundaries in L0
    yrange = [-1.25 1.25];  % y boundaries in L0

    %% in this test case, the condition number is almost
    % twice as large, but the error is somehow less?
    % xrange = [-1 1];  % x boundaries in L0
    % yrange = [-1 1];  % y boundaries in L0

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
    thresh = 1;
    %[L,U,P,Q] = lu(A);
     %%more stable factorization (since we're using it repeatedly in born)
    [L,U,P,Q,D] = lu(A, thresh);
    toc

    tic
    Ez_for = reshape(P.'*(U\(L\(Q.'*(D\b)))),N);
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
    dEz_for{1}=Ez_for; % dEz_for... the derivative of dEz...usualy computed by adjoint.

    %II1b, II2b is the point where we are optimizing
    % so are we calculating the series only for one point?
    % yes, because we are only optimizing one point, so we only 
    % need to calculate the field at that point for a given deltaEpsilon

    S(1)=dEz_for{1}(II1b,II2b);

    %% SWITCH BETWEEN N_order = 20
    %% N_order = 7 should produce convergent...
    N_order=13;

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
        %new_field = reshape((P.'*(U\(L\(Q.'* bprime)))),N);
        new_field = reshape(P.'*(U\(L\(Q.'* (D\bprime)))),N);
        dEz_for{iorder} = new_field;
        %[dEz_for{iorder}] = solveTM(wvlen, xrange, yrange, eps_str_old, Mz, Npml);

    end;
    toc

    figure();
    imagesc(real(Mz_LS_cell{2})); colorbar;

    % I suspect this is the search grid, but not sure
    Shanks_cell{1} = dEz_for{1};

    alphav = [0.9]
    ialpha = 1

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

    for iorder=1:N_order-1
        e{2}(iorder)=1/(S(iorder+1)-S(iorder));
        eShanks{2}{iorder} = 1./(Shanks_cell{iorder+1} - Shanks_cell{iorder});
        %eShanks{2}{iorder}(zone==0) = 0;
    end

    for ip=3:N_order
        eShanks{ip} = {};
        for iorder=1:N_order-ip+1
            % this looks much like the shanks transformation
            e{ip}(iorder)=e{ip-2}(iorder+1)+1/(e{ip-1}(iorder+1)-e{ip-1}(iorder)); 
            eShanks{ip}{iorder} = eShanks{ip-2}{iorder+1} + 1./(eShanks{ip-1}{iorder+1} - eShanks{ip-1}{iorder});
            %eShanks{ip}{iorder}(zone == 0) = 0;
        end

    end



    %% validate prediction
    new_eps = deeps_str+eps_str_old;
    Ez_validation = solveTM(wvlen, xrange, yrange, new_eps, Mz0, Npml);

    figure();
    subplot(1,3,1)
    imagesc(real(Ez_validation)); colorbar
    valcaxis = caxis;
    subplot(1,3,2)
    imagesc(real(eShanks{N_order}{1})); colorbar;
    caxis(valcaxis)
    subplot(1,3,3)
    imagesc(abs(eShanks{N_order}{1}-Ez_validation)); colorbar;
    title('difference')

    % shanks_order = N_order-6;

    for shanks_order = 3:2:N_order
    difference = abs(eShanks{shanks_order}{1}-Ez_validation);

    bounds = [0.0]; %[0, 1e-3, 0.01, 0.02, 0.03]
    colors = ['r', 'g', 'b', 'c', 'm'];
    figure();
    for b = 1:length(bounds)
        difference_mask = difference> bounds(b);
        [row,col] = find(difference_mask == 1);


        %% analyze the shanks series for these special points

        pointDiffs = [];
        fields = [];
        for i = 1:length(row)

           eShanks{shanks_order-1}{1}(row(i), col(i));
           eShanks{shanks_order-1}{2}(row(i), col(i));

           f1 = eShanks{shanks_order-2}{1+1}(row(i), col(i));
           point_diff = eShanks{shanks_order-1}{2}(row(i), col(i)) - eShanks{shanks_order-1}{1}(row(i), col(i));

           val = f1+1/point_diff;
           actual_shanks = eShanks{shanks_order}{1}(row(i), col(i));
           actual = Ez_validation(row(i), col(i));
           pdifference = abs(actual-actual_shanks);
           pointDiffs(i,1) = point_diff;
           pointDiffs(i,2) = pdifference; 
           fields(i) = actual;
        end

        %% this is the most illustrative figure...
        % there appears to be an exponential increase in the error 
        % based on the magnitude of pointDiff
        subplot(1,2,1)
        plot(abs(pointDiffs(:,1)), pointDiffs(:,2), '.', 'MarkerSize',5)
        xlabel('denominator difference')
        ylabel('shanks error against actual field')
        hold on;
        %% we should compare it to the general population of errors as well
        xlim([0,100])

        subplot(1,2,2);
        if(b==1)
            imagesc(difference);
            hold on;
        end
        plot(col, row, '.', 'MarkerFaceColor', colors(b))
        hold on;
    end
    drawnow();
    end
    data(wi,1) = condest(A);
    data(wi,2) = median(pointDiffs(:,2));

end

figure(); plot(data(:,1), data(:,2), '.', 'MarkerSize',20);
xlabel('condition of A')
ylabel('mean error')