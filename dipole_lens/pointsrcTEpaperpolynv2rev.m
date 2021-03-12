clear all;
close all;
clc;

%% SCRIPT IS TO A BRUTE FORCE or RIGOROUS ASSESSMENT OF LEARNING RATE

%%
spparms('spumoni',0)

figure(1)
%figure(2)
figure(3)
figure(4)
drawnow
%pause(25)

eps_0 = 8.85*10^-12;
mu_0 = 4*pi*10^-7; 
Z0=sqrt(mu_0/eps_0);

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 1.55;%1.565;  % wavelength in L0
xrange0 = [-2-0.025 2+0.025];  % x boundaries in L0
yrange0 = [-2-0.025 2+0.025];  % y boundaries in L0
N0 = 1.0.*[81 81];  % [Nx Ny]
Npml = [30 30];  % [Nx_pml Ny_pml]
Nx = N0(1); Ny = N0(2);
k0=(2*pi)/wvlen;

[YY0,XX0]=meshgrid(linspace(yrange0(1)+0.025,yrange0(2)-0.025,N0(2)),linspace(xrange0(1)-0.025,xrange0(2)-0.025,N0(1)));

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange0, yrange0, N0, Npml);  % domain is expanded to include PML

[YY,XX]=meshgrid(linspace(yrange(1)+0.025,yrange(2)-0.025,N(2)),linspace(xrange(1)+0.025,xrange(2)-0.025,N(1)));

TT0b=(XX<(xrange(2)-Lpml(1))).*(XX>(xrange(1)+Lpml(1))).*(YY<(yrange(2)-Lpml(2))).*(YY>(yrange(1)+Lpml(2)));
%TTopt=(XX<(xrange(2)-Lpml(1))).*(XX>(xrange(1)+Lpml(1))).*(YY<(yrange(2)-Lpml(2))).*(YY>(yrange(1)+Lpml(2)));
TTopt=(XX<=1.55).*(XX>=-1.55).*(YY<=1.55).*(YY>=-1.55);


XX0b=XX(TT0b==1);
YY0b=YY(TT0b==1);

XXopt=XX(TTopt==1);
YYopt=YY(TTopt==1);

k0=2*pi/wvlen;


%% Set up the permittivity.
ff=0.42;
GGG=1;%0.25;%0.015;%0.25;%0.01;
epsilon = 1;
numCells = 10;
cellsize = Nx/numCells;
%load('struct.txt');
%eps_str=(1.5)^2.*ones(N);
%eps_str(21:83,21:83)=eps_str(21:83,21:83)+struct.'.*(2.5^2-1);
%ky=-25*pi/(yrange0(2)-yrange0(1))+(1:51).*pi/(yrange0(2)-yrange0(1));


irandomfini=0;
eps_str = ((1.5)^2+1.1*GGG*0).*ones(N);
[kx,ky]=meshgrid(linspace(-10*pi/(xrange0(2)-xrange0(1)),10*pi/(xrange0(2)-xrange0(1)),51),linspace(-10*pi/(yrange0(2)-yrange0(1)),10*pi/(yrange0(2)-yrange0(1)),51));
randxy=exp(i.*(rand(51,51).*2*pi));
    
    
FF=0;
for ikx=1:51,
    for iky=1:51,
FF=FF+randxy(ikx,iky).*exp(i.*kx(ikx,iky).*XX+1i.*ky(ikx,iky).*YY);
end;
end;
FF=real(FF);
FF=((2.5)^2-(1.5)^2).*(FF-min(min(FF)))./(max(max(FF))-min(min(FF)));

%FF(FF>=((2.5)^2-(1.5)^2-1.1*GGG)/2)=((2.5)^2-(1.5)^2-1.1*GGG);
%FF(FF<((2.5)^2-(1.5)^2-1.1*GGG)/2)=0;

%FF=(2.5^2-1.5^2).*round(rand(size(XX,1),size(XX,2)))

eps_str(TTopt==1)=eps_str(TTopt==1)+FF(TTopt==1);%3.*(rand(size(eps_str(TTopt==1))));


[C1a,I1a]=min(sqrt((XX-0).^2+(YY-0).^2));
[C2a,II2a]=min(C1a);
II1a=I1a(II2a);

[C1b,I1b]=min(sqrt((XX-0).^2+(YY-1.55).^2));
[C2b,II2b]=min(C1b);
II1b=I1b(II2b);

[C1c,I1c]=min(sqrt((XX+3).^2+(YY+4).^2));
[C2c,II2c]=min(C1c);
II1c=I1c(II2c);

Mz = zeros(N);
Mz(II1a,II2a)=1;

%EPSTOT=2*ones(N);
%Mz(YY==min(min(YY0b)))=exp(-0.5.*(XX(YY==min(min(YY0b)))+0).^2);

[Ez_m, Hx_m, Hy_m, A, omega,b] = solveTM(L0,wvlen, xrange, yrange, eps_str, Mz, Npml);

GGauto=zeros(N);
iiter=1;
alphafix=1;
%binaire=1;
%newbinaire=1;
binaire=0;cptbinaire=0;newbinaire=0;
cptadapt=1;
gradient_history = cell(1);
for iiter=1:3000,%7,%200;50;%2000,%75;%2000,%while irandomfini==0,%for iiter=1:2200,

     tic
    for iadj=1:2,
        %% Set up the magnetic current source density.
        Mz = zeros(N);
        %ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
        %Mz(ind_src(1), ind_src(2)) = 1;

        %if iadj==1,Mz(YY==max(max(YY0b)))=exp(-1.*(XX(YY==max(max(YY0b)))+1).^2);else,Mz(YY==min(min(YY0b)))=exp(-1.*(XX(YY==min(min(YY0b)))-1).^2);end;
        if iadj==1,Mz(II1a,II2a)=1./(-i*k0*Z0*0.05^2);
        else,
            Mz(II1b,II2b)=1./(-i*k0*Z0*0.05^2);
            %Mz(YY==min(min(YY0b)))=exp(-0.5.*(XX(YY==min(min(YY0b)))+0).^2);
    end;




    if iadj==1,


        disp('resol forward')    
        [Ez_for, Hx_for, Hy_for, A, omega,b] = solveTM(L0,wvlen, xrange, yrange, eps_str, Mz, Npml);

        else,
        tic
        disp('resol adj')
        [Ez_adj, Hx_adj, Hy_adj, A, omega,b] = solveTM(L0,wvlen, xrange, yrange, eps_str, Mz, Npml);

    end;

    end;


    Ez_for_old=Ez_for;

    Objective(iiter)=log10((abs(Ez_for(II1b,II2b)))^2)+2.369+0.04936;%./imag(Ez_for(II1a,II2a));
    Objectiveadapt(cptadapt)=log10((abs(Ez_for(II1b,II2b)))^2)+2.369+0.04936;%


    Timecpt(iiter)=toc;

    figure(3),
    plot(Objective)

    dFOM=real((2.5^2-1.5^2).*Ez_for.*conj(Ez_for(II1b,II2b)).*Ez_adj.*exp(-i.*0*angle(Ez_m(II1a,II2a))));
    dFOM(dFOM>0)=dFOM(dFOM>0).*(2.5^2-eps_str(dFOM>0)).^(2)./((2.5^2-1.5^2).^(2));
    dFOM(dFOM<0)=dFOM(dFOM<0).*(eps_str(dFOM<0)-1.5^2).^(2)./((2.5^2-1.5^2).^(2));

    if iiter>1, eps_str_old2=eps_str_old;end
    eps_str_old=eps_str;


    facmin=min(min((2.5^2-eps_str((TTopt.*(dFOM>0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1))./dFOM((TTopt.*(dFOM>0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)),min((eps_str((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)-1.5^2)./(-dFOM((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1))));
   
    eps_strnorm=eps_str;
    eps_strnorm(TTopt==1)=eps_strnorm(TTopt==1)+dFOM(TTopt==1).*(GGG).*facmin;%0.005

    eps_str((TTopt.*(eps_str>2.5^2)==1))=2.5^2;
    eps_str((TTopt.*(eps_str<1.5^2)==1))=1.5^2;

    deeps_str=eps_strnorm-eps_str_old;
    gradient_history{iiter} = deeps_str;
    Ez_proj{iiter}=Ez_for;
    dEz_for{1}=Ez_for;
    S(1)=dEz_for{1}(II1b,II2b);
    
    %% something's wrong with the N_order...here, it has to be way too large?
    N_order=3;%5
    for iorder=2:N_order,
        Mz(TTopt==1)=1./(-i*k0*Z0)*k0^2*deeps_str(TTopt==1).*dEz_for{iorder-1}(TTopt==1);
       [dEz_for{iorder}, Hx_forin, Hy_forin, Ain, omegain,bin] = solveTM(L0,wvlen, xrange, yrange, eps_str_old, Mz, Npml);
    end;

  
    alphav=[1e-6:1e-6:9e-6 1e-5:1e-5:9e-5 1e-4:1e-4:9e-4 1e-3:1e-3:9e-3 1e-2:1e-2:9e-2 1e-1:1e-1:1];%0:0.0001:1;

    for ialpha=1:numel(alphav),
% 
%         for iorder=2:N_order,
%             Ez_proj{iiter}=Ez_proj{iiter}+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder};
%             S(iorder)=S(iorder-1)+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder}(II1b,II2b);
%         end;
% 
%         e{1}=S;
%         for iorder=1:N_order-1,
%             e{2}(iorder)=1/(S(iorder+1)-S(iorder));
%         end;
%         for ip=3:N_order,
%             for iorder=1:N_order-ip+1,
%                 e{ip}(iorder)=e{ip-2}(iorder+1)+1/(e{ip-1}(iorder+1)-e{ip-1}(iorder));    
%             end;
%         end;
%         prediction(ialpha)=log10((abs(S(N_order))).^2)+2.369+0.04936;
%         
        %% CHECK THIS
        eps_str_brute=eps_str_old;
        eps_str_brute(TTopt==1)=eps_str_brute(TTopt==1)+dFOM(TTopt==1).*(alphav(ialpha)).*facmin;%0.005
        Mz = zeros(N);
        Mz(II1a,II2a)=1./(-i*k0*Z0*0.05^2);
        [Ez_for_brute, Hx_for_brute, Hy_for_brute, Ain, omegain,bin] = solveTM(L0,wvlen, xrange, yrange, eps_str_brute, Mz, Npml);
        obj_brute(ialpha)=log10((abs(Ez_for_brute(II1b,II2b)))^2)+2.369+0.04936;

    end;

    UUUint=1e-6:1e-6:1;
    %VVVint=interp1(alphav,full(prediction),UUUint,'spline');
    %WWW(iiter)=UUUint(find(VVVint==max(VVVint)));
    VVVbrute = interp1(alphav,full(obj_brute),UUUint,'spline');
    WWWbrute(iiter)=UUUint(find(VVVbrute==max(VVVbrute)));
    
    %% use the brute force version
    alphamax(iiter)=WWWbrute(iiter);

    eps_str(TTopt==1)=eps_str(TTopt==1)+dFOM(TTopt==1).*(alphamax(iiter)).*facmin;%0.005
    predictedObjective=log10((abs(Ez_proj{iiter}(II1b,II2b)))^2)+2.369+0.04936;

  
    drawnow
    iiter=iiter+1;
    iiter
    cptadapt=cptadapt+1;

    
end;

save("alphamax_data.mat", "alphamax", "gradient_history","Objective");



