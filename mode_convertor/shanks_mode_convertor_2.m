clear all;
close all;
clc;

%%
%spparms('spumoni',0)

%load('epstest.txt')

load('modeTE0re.txt')
load('modeTE0im.txt')

load('modeTE1re.txt')
load('modeTE1im.txt')

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

XX0b=XX(all([XX<(xrange(2)-Lpml(1)) XX>(xrange(1)+Lpml(1)) YY<(yrange(2)-Lpml(2)) YY>(yrange(1)+Lpml(2))]));
YY0b=YY(all([XX<(xrange(2)-Lpml(1)) XX>(xrange(1)+Lpml(1)) YY<(yrange(2)-Lpml(2)) YY>(yrange(1)+Lpml(2))]));

TT0b=(XX<(xrange(2)-Lpml(1))).*(XX>(xrange(1)+Lpml(1))).*(YY<(yrange(2)-Lpml(2))).*(YY>(yrange(1)+Lpml(2)));
%TTopt=(XX<(xrange(2)-Lpml(1))).*(XX>(xrange(1)+Lpml(1))).*(YY<(yrange(2)-Lpml(2))).*(YY>(yrange(1)+Lpml(2)));
%TTopt=(XX<=1.55).*(XX>=-1.55).*(YY<=1.55).*(YY>=-1.55);
TTopt=(XX<=1).*(XX>=-1).*(YY<=1).*(YY>=-1);

TTmode1=(XX>=-1.55-0.025).*(XX<=1.55+0.025).*(YY<=-1.49).*(YY>=-1.51);%(YY==(min((YY(TTopt==1)))-0.15));
TTmode2=(XX>=-1.55-0.025).*(XX<=1.55+0.025).*(YY<=1.51).*(YY>=1.49);%(YY==(max((YY(TTopt==1)))+0.15));


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


for irandom=1:1,%10,
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
FF=((2.5)^2-(1.5)^2-1.1*GGG*0).*(FF-min(min(FF)))./(max(max(FF))-min(min(FF)));

%FF(FF>=((2.5)^2-(1.5)^2-1.1*GGG*0)/2)=((2.5)^2-(1.5)^2-1.1*GGG*0);
%FF(FF<((2.5)^2-(1.5)^2-1.1*GGG*0)/2)=0;

%FF=(2.5^2-1.5^2).*round(rand(size(XX,1),size(XX,2)))

%eps_str(TTopt==1)=eps_str(TTopt==1)+3.*(rand(size(eps_str(TTopt==1))));



%eps_str(TTopt==1)=eps_str(TTopt==1)+FF(TTopt==1);

% load('soluce.txt');
% 
% for ix=1:size(XX,1),
%     for iy=1:size(XX,2),
%         
%         if any((abs(soluce(:,1)-XX(ix,iy))+abs(soluce(:,2)-YY(ix,iy)))<0.01),
%             
%             eps_str(ix,iy)=2.5.^2;
%         end;
%         
%         %if all([YY(ix,iy)>1.55 YY(ix,iy)<max(max(YY0b))+1 abs(XX(ix,iy))<0.15]),
%         %    eps_str(ix,iy)=2.5.^2;
%         %end;
%         
%     end;
% end;

eps_str(((XX<=0.5).*(XX>=-0.5))==1)=2.5^2;
%eps_str(((XX<=0.5).*(XX>=-0.5).*(YY<-1))==1)=2.5^2;
%eps_str(((XX<=0.5).*(XX>=-0.5).*(YY>=1))==1)=2.5^2;

%eps_str(TTopt==1)=2^2;%+0.001;%eps_str(TTopt==1)+FF(TTopt==1);%1.5^2+0.001;%



%eps_str=epstest;


%savefig('papersrclinbinini.fig')


%eps_str(TTopt.*(XX>-(2.*YY+15)./5).*(XX<(5-2.*YY)./5)==1)=eps_str(TTopt.*(XX>-(2.*YY+15)./5).*(XX<(5-2.*YY)./5)==1)+3;

%eps_str(TTopt==1)=eps_str(TTopt==1)+1.5;

figure(15)
pcolor(XX,YY,eps_str);shading('interp'),colorbar,axis([min(min(XX0b)) max(max(XX0b)) min(min(YY0b)) max(max(YY0b))]);title('eps_str')
drawnow

Mz = zeros(N);
Mz(TTmode1==1)=(modeTE0re+1i.*modeTE0im);
[Ez_inc, Hx_inc, Hy_inc, Atest] = solveTM(L0,wvlen, xrange, yrange, eps_str, Mz, Npml);
Pinc=0.5*sum(real((Ez_inc(TTmode2==1).*conj(Hx_inc(TTmode2==1)))))*0.05;

figure(1)
pcolor(XX,YY,((real(Ez_inc))));title('real(Einc)');shading('interp'),colorbar,axis([-1.55 1.55 -1.55 1.55]),xlabel('X(um)'),ylabel('Z(um)')
%caxis([-1 4])
drawnow, pause(5)


Mz = zeros(N);
Mz(TTmode1==1)=(modeTE1re+i.*modeTE1im);
[Ez_mode, Hx_mode, Hy_moden] = solveTM(L0, wvlen, xrange, yrange, eps_str, Mz, Npml);
Pmode=0.5*sum(real((Ez_mode(TTmode2==1).*conj(Hx_mode(TTmode2==1)))))*0.05;

figure(1)
pcolor(XX,YY,((real(Ez_mode))));title('real(Emode)');shading('interp'),colorbar,axis([-1.55 1.55 -1.55 1.55]),xlabel('X(um)'),ylabel('Z(um)')

%eps_str(TTopt==1)=1.50001^2;
eps_str(TTopt==1)=2^2;

%eps_str(TTopt==1)=eps_str(TTopt==1)+rand(1681,1).*3.5;


alphafix=1;
GGauto=zeros(N);
iiter=1;
%binaire=1;
%newbinaire=1;
binaire=0;cptbinaire=0;newbinaire=0;
IITERS = 550;
for iiter=1:IITERS,%7,%200;50;%2000,%75;%2000,%while irandomfini==0,%for iiter=1:2200,

    %if iiter>3,if Objective(iiter-1)<Objective(iiter-2),GGG=0.5*GGG;end;end;
    
    %if iiter==9987,return,end;
    tic
for iadj=1:2,
%% Set up the magnetic current source density.
Mz = zeros(N);
%ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
%Mz(ind_src(1), ind_src(2)) = 1;

%if iadj==1,Mz(YY==max(max(YY0b)))=exp(-1.*(XX(YY==max(max(YY0b)))+1).^2);else,Mz(YY==min(min(YY0b)))=exp(-1.*(XX(YY==min(min(YY0b)))-1).^2);end;
if iadj==1,Mz(TTmode1==1)=(modeTE0re+i.*modeTE0im);
else,
    Mz(TTmode2==1)=(modeTE1re+i.*modeTE1im);
    %Mz(YY==min(min(YY0b)))=exp(-0.5.*(XX(YY==min(min(YY0b)))+0).^2);
end;


%abs(YY-4)<0.05)=1;%exp(-1.*(abs(XX(abs(YY-4)<0.05))).^2);%exp(-i.*k0.*YY(and(abs(XX)<3,abs(YY-4)<0.01)));

%Mz=exp(i.*k0.*YY);


%% Solve TE equations.

% tic
% [Hz, Ex, Ey,A,omega,b, Sxf, Dxf, Dyf, sxf, syf,trun] = solveTE(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
% toc

if iadj==1,

disp('resol forward')    
% [Ez_for_check, Hx_for_check, Hy_for_check] = solveTM(L0, wvlen, xrange, yrange, eps_str, Mz, Npml);
 [Ez_for, Hx_for, Hy_for, A,omega, b] = ...
 solveTM_dirichlet(L0, wvlen, xrange, yrange, eps_str, Mz, Npml);
 
else,
tic
disp('resol adj')
%[Ez_adj, Hx_adj, Hy_adj] = solveTM(L0, wvlen, xrange, yrange, eps_str, Mz, Npml);
[Ez_adj, Hx_adj, Hy_adj] = solveTM_dirichlet(L0, wvlen, xrange, yrange, eps_str, Mz, Npml);

end;
    
end;

figure(15)
pcolor(XX,YY,eps_str);shading('interp'),colorbar,axis([min(min(XX0b)) max(max(XX0b)) min(min(YY0b)) max(max(YY0b))]);title('eps_str')


toc
figure(1)
%pcolor(XX,YY,log10((abs(Ez_for)).^2)+2.369);title('log_{10}(E^2)');shading('interp'),colorbar,axis([-1.55 1.55 -1.55 1.55].),caxis([-1 4]),xlabel('X(um)'),ylabel('Z(um)')
pcolor(XX,YY,((real(Ez_for)./max(max(real(Ez_inc))))));title('real(E)');shading('interp'),colorbar,axis([-1.55 1.55 -1.55 1.55]),xlabel('X(um)'),ylabel('Z(um)')
frame = getframe(figure(1));
im1{iiter} = frame2im(frame);


if iiter>1,
figure(10)
pcolor(XX,YY,((real(Ez_for_old))));title('real(E_{old})');shading('interp'),colorbar,axis([-1.55 1.55 -1.55 1.55]),caxis([-1 4]),xlabel('X(um)'),ylabel('Z(um)')
    
end;


1/(4*sqrt(Pmode*Pinc))*(0.05*(sum(Ez_for(TTmode2==1).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*(Hx_for(TTmode2==1))))),

if iiter>1,

figure(50),
plot(XX(TTmode2==1),real(SEz{1,N_order}),'r');title(sprintf('Ez projection %g',iiter)),hold on
plot(XX(TTmode2==1),real(Ez_for(TTmode2==1)),'k');
plot(XX(TTmode2==1),real(Ez_for_old(TTmode2==1)),'k--');hold off,xlabel('X(µm)')

figure(51),
plot(XX(TTmode2==1),real(eHx{1,N_order}),'r');title(sprintf('Hx projection %g',iiter)),hold on
plot(XX(TTmode2==1),real(Hx_for(TTmode2==1)),'k');
plot(XX(TTmode2==1),real(Hx_for_old(TTmode2==1)),'k--');hold off,xlabel('X(µm)')

figure(52),
plot(XX(TTmode2==1),real(eHy{1,N_order}),'r');title(sprintf('Hy projection %g',iiter)),hold on
plot(XX(TTmode2==1),real(Hy_for(TTmode2==1)),'k');
plot(XX(TTmode2==1),real(Hy_for_old(TTmode2==1)),'k--');hold off,xlabel('X(µm)')
%figure(60)
%pcolor(XX,YY,log10(abs(Ez_proj{iiter-1}-Ez_for)./abs(Ez_for)));title('relative error Ez');shading('interp'),colorbar,axis([-1.55 1.55 -1.55 1.55].),xlabel('X(um)'),ylabel('Z(um)')
if 0,%iiter==2,
    %1/(4*sqrt(Pmode*Pinc))*(0.05*(sum(Ez_for(TTmode2==1).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*(Hx_for(TTmode2==1))))),
    %1/(16*(Pmode*Pinc))*(0.05*abs(sum(Ez_for(TTmode2==1).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*(Hx_for(TTmode2==1)))))^2,
    return,pause(10),end;    
end;

Ez_for_old=Ez_for;
Hx_for_old=Hx_for;
Hy_for_old=Hy_for;

%Objective(iiter)=(abs(sum(Ez_for(YY==min(min(YY0b))).*conj(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2)))))^2./((sum((abs(Ez_for(YY==min(min(YY0b))))).^2)).*(sum((abs(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2))).^2)));


%Objective(iiter)=log10((abs(Ez_for(II1b,II2b)))^2)+2.369;%./imag(Ez_for(II1a,II2a));

Objective(iiter)=1/(16*Pmode*Pinc)*(0.05*abs(sum(Ez_for(TTmode2==1).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*(Hx_for(TTmode2==1))))).^2;
Aconj(iiter)=1/(4*Pmode*Pinc)*(0.05*(sum(Ez_for(TTmode2==1).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*(Hx_for(TTmode2==1)))));

%(yrange0(2)-yrange0(1))/Ny.*

Timecpt(iiter)=toc;

if iiter>1,
   if Objective(iiter)<Objective(iiter-1),
alphafix=alphafix*0.9;       
   end;
end;

%A=real(Ez_for.*conj(Ez_for(II1b,II2b)).*Ez_adj).*exp(-i.*angle(Ez_m(II1a,II2a)));
%A=sum(-Ez_for(XX==max(max(XX0b))).*conj(Hy_m(XX==max(max(XX0b))))-conj(-Ez_m(XX==max(max(XX0b)))).*(Hy_for(XX==max(max(XX0b))))+Ez_for(YY==min(min(YY0b))).*conj(Hx_m(YY==min(min(YY0b))))+conj(-Ez_m(YY==min(min(YY0b)))).*(Hx_for(YY==min(min(YY0b)))));

%B=Ez_m(II1a))));

%dFOM=(eps_str-1).*real(Ez_for.*conj(sum(Ez_for(YY==min(min(YY0b))).*conj(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2)))).*Ez_adj.*exp(-i.*angle(sum(Ez_adj(YY==min(min(YY0b))).*conj(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2)))))./((sum((abs(Ez_for(YY==min(min(YY0b))))).^2)).*(sum((abs(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2))).^2))));

%dFOM=real(Ez_for.*conj(sum(Ez_for(YY==min(min(YY0b))).*conj(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2)))).*Ez_adj.*exp(-i.*angle(sum(Ez_adj(YY==min(min(YY0b))).*conj(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2)))))./((sum((abs(Ez_for(YY==min(min(YY0b))))).^2)).*(sum((abs(exp(-1.*(XX(YY==min(min(YY0b)))-1).^2))).^2))));

%dFOM=real(Ez_for.*conj(Ez_for(II1b,II2b)./imag(Ez_for(II1a,II2a))).*Ez_adj.*exp(-i.*angle(Ez_m(II1a,II1b)./imag(Ez_m(II1a,II2a)))));%.*exp(-i.*angle(B)));

figure(3),
plot(log10(1-Objective))

% dFOMp=real((2.5^2-1.5^2).*Ez_for./(1+GGauto.*(2.5^2-1.5^2)).*conj(Ez_for(II1b,II2b)).*Ez_adj.*exp(-i.*angle(Ez_m(II1a,II2a))));
% dFOMn=real(-(2.5^2-1.5^2).*Ez_for./(1-GGauto.*(2.5^2-1.5^2)).*conj(Ez_for(II1b,II2b)).*Ez_adj.*exp(-i.*angle(Ez_m(II1a,II2a))));
% dFOM=(dFOMp>=dFOMn).*(dFOMp>=0).*dFOMp-(dFOMp<dFOMn).*(dFOMn>=0).*dFOMn;

%dFOM=real((2.5^2-1.5^2).*Ez_for.*conj(Ez_for(II1b,II2b)).*Ez_adj.*exp(-i.*0*angle(Ez_m(II1a,II2a))));
dFOM=real(Ez_for.*conj(Aconj(iiter)).*Ez_adj.*exp(i.*21.13*pi/180));

dFOM(dFOM>0)=dFOM(dFOM>0).*(2.5^2-eps_str(dFOM>0)).^(1)./((2.5^2-1.5^2).^(1));
dFOM(dFOM<0)=dFOM(dFOM<0).*(eps_str(dFOM<0)-1.5^2).^(1)./((2.5^2-1.5^2).^(1));

%dFOM(TTopt==0)=0;
%[C1b,I1b]=max(abs(dFOM));
%[C2b,II2b]=max(C1b);
%II1b=I1b(II2b);

if ~binaire,%iiter<=70,

eps_str((TTopt.*(eps_str>2.5^2)==1))=2.5^2;
eps_str((TTopt.*(eps_str<1.5^2)==1))=1.5^2;

    eps_str_old=eps_str;

facmin=min(min((2.5^2-eps_str((TTopt.*(dFOM>0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1))./dFOM((TTopt.*(dFOM>0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)),min((eps_str((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)-1.5^2)./(-dFOM((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1))));
%facmin=min(min((2.5^2-1.5^2)./dFOM((TTopt.*(dFOM>0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)),min((1.5^2-2.5^2)./dFOM((TTopt.*(dFOM<0).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1)));

%if iiter==6,GGG=1;end;

%eps_str(TTopt==1)=eps_str(TTopt==1)+dFOM(TTopt==1).*(GGG+0.*iiter/1000)./max(max(abs(dFOM(((TTopt).*(eps_str>(1.5^2+GGG)).*(eps_str<(2.5^2-GGG)))==1))));%0.005
%eps_str(TTopt==1)=eps_str(TTopt==1)+dFOM(TTopt==1).*(GGG+0.*iiter/1000)./max(max(abs(dFOM(((TTopt).*(eps_str>(1.5^2)).*(eps_str<(2.5^2)))==1))));%0.005

%remettre celui la
eps_strnorm=eps_str;
eps_strnorm(TTopt==1)=eps_strnorm(TTopt==1)+dFOM(TTopt==1).*(GGG).*facmin;%0.005

figure(16)
pcolor(XX,YY,eps_str);shading('interp'),colorbar,title('eps new'),axis([-1.55 1.55 -1.55 1.55]);

figure(17)
pcolor(XX,YY,eps_str-eps_str_old);shading('interp'),title('deeps'),colorbar,axis([-1.55 1.55 -1.55 1.55]);


deeps_str=eps_strnorm-eps_str_old;
if iiter==49
    figure(40),pcolor(XX,YY,deeps_str),axis([-1.55 1.55 -1.55 1.55]),shading('interp'),colorbar,title('\Delta eps'),xlabel('X(um)'),ylabel('Z(um)')
    figure(41),surf(XX,YY,deeps_str-(eps_str-eps_str_old)),axis([-1.55 1.55 -1.55 1.55]),shading('interp'),colorbar,title('\Erreur Delta eps'),xlabel('X(um)'),ylabel('Z(um)')
end;

clear eEz eHx eHy SEz SHx SHy
Ez_proj{iiter}=Ez_for;

dEz_for{1}=Ez_for;
dHx_for{1}=Hx_for;
dHy_for{1}=Hy_for;

S(1)=Aconj(iiter);%dEz_for{1}(II1b,II2b);
SEz{1}=dEz_for{1}(TTmode2==1);
SHx{1}=dHx_for{1}(TTmode2==1);
SHy{1}=dHy_for{1}(TTmode2==1);
% SEz{1}=(SEz{1}).';
% SHx{1}=(SHx{1}).';
% SHy{1}=(SHy{1}).';

%% what should this be properl
N_order=3;%5 1.3475e-05 + 9.0151e-07i 0.9109 + 0.4477i 1.0163 + 0.2631i 0.5972 + 0.0144i 0.5262 + 0.0984i 0.5503 + 0.1912i 0.6029 + 0.2269i 0.6858 + 0.1489i 0.6349 + 0.0853i 0.5670 + 0.1367i 0.5929 + 0.1805i 0.6295 + 0.1680i 0.6316 + 0.1446i 0.6132 + 0.1335i 0.5987 + 0.1455i 0.6074 + 0.1612i 0.6208 + 0.1565i 0.6189 + 0.1453i 0.6108 + 0.1448i 0.6085 + 0.1505i 0.6122 + 0.1538i
%1.3475e-05 + 9.0151e-07i 0.7316 + 0.3376i 0.6884 + 0.1624i 0.6051 + 0.1823i 0.6203 + 0.1461i 0.6125 + 0.1496i 0.6131 + 0.1501i 0.6130 + 0.1501i 0.6130 + 0.1501i 0.6130 + 0.1501i 0.6130 + 0.1501i 
%1 1 0.9470 0.9740 0.9260 0.9300 0.9300 0.9300 0.9300 0.9300
%1 1 0.8490 0.7760 0.8260 1 1 1 0.8230 0.8830 1 1 0.9050 0.8700 0.9110 1 1 0.9040 0.8990 0.9320 1 1 0.9130 0.9150 0.9410 1 0.9330 0.9210 0.9240 0.9360 0.9390 0.9290 0.9250 0.9280 0.9330 0.9320 0.9290 0.9280 0.9290 0.9310    

%1 0.5060 0.5040 0.5030 0.5030 0.5030 0.5030 0.5030 0.5030 0.5030 
%1 1 0.4680 0.46prediction(ialpha)=1/(16*Pmode*Pinc)*(0.05*abs(sum((eEz{1,N_order}).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((eHx{1,N_order}))))).^2;30 0.4950 0.5120 0.5080 0.5030 0.5030 0.5030 0.5040 0.5040 0.5030 0.5030 0.5030 0.5030 0.5030 0.5030 0.5030 0.5030   
for iorder=2:N_order,
    Mz = zeros(N);
    Mz(TTopt==1)=1./(-i*k0*Z0)*k0^2*deeps_str(TTopt==1).*dEz_for{iorder-1}(TTopt==1);
   [dEz_for{iorder}, dHx_for{iorder}, dHy_for{iorder}] = solveTM(L0,wvlen, xrange, yrange, eps_str_old, Mz, Npml);
end;


%     for ix=1:size(Ez_for,1),
%     for iy=1:size(Ez_for,2),
%         
% dEz_for{iorder}(ix,iy)=sum(deeps_str(TTopt==1)*k0^2.*0.05^2*i/4.*besselh(0,k0*1.5*sqrt((XX(ix,iy)-XX(TTopt==1)).^2+(YY(ix,iy)-YY(TTopt==1)).^2)).*dEz_for{iorder-1}(ix,iy));    
% 
%     end;
%     end;

%alphav=[1e-6 1e-5 1e-4 0.001:0.001:1];%0:0.001:1;
alphav=[1e-6 1e-5 1e-4 1e-3 0.01:0.01:1];
for ialpha=1:numel(alphav),
    
for iorder=2:N_order,
    Ez_proj{iiter}=Ez_proj{iiter}+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder};

    Aconjorder=1/(4*Pmode*Pinc)*(0.05*(sum(dEz_for{iorder}(TTmode2==1).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*(dHx_for{iorder}(TTmode2==1)))));
    S(iorder)=S(iorder-1)+(alphav(ialpha)).^(iorder-1)*Aconjorder;%dEz_for{iorder}(II1b,II2b);

    SEz{iorder}=SEz{iorder-1}+(alphav(ialpha)).^(iorder-1)*dEz_for{iorder}(TTmode2==1);%dEz_for{1}(TTmode2==1);
    SHx{iorder}=SHx{iorder-1}+(alphav(ialpha)).^(iorder-1)*dHx_for{iorder}(TTmode2==1);
    SHy{iorder}=SHy{iorder-1}+(alphav(ialpha)).^(iorder-1)*dHy_for{iorder}(TTmode2==1);

    %figure(90+(iiter-1)*20+iorder-2)
%pcolor(XX,YY,((real(dEz_for{iorder}))));title(sprintf('log_{10}(E_{extrap}^2) order=%g)',iorder-1));shading('interp'),colorbar,axis([-1.55 1.55 -1.55 1.55].),caxis([min(min(real(dEz_for{iorder}))) max(max(real(dEz_for{iorder})))]),xlabel('X(um)'),ylabel('Z(um)')

end;


for iorder=1:N_order,
    SEz{iorder}=full(SEz{iorder});
    SHx{iorder}=full(SHx{iorder});
    SHy{iorder}=full(SHy{iorder});
end;

e{1}=S;

for iorder=1:N_order,%ix=1:sum(sum(TTmode2==1)),

    eEz{iorder,1}=SEz{iorder};
    eHx{iorder,1}=SHx{iorder};
    eHy{iorder,1}=SHy{iorder};
end;

for iorder=1:N_order-1,
    e{2}(iorder)=1/(S(iorder+1)-S(iorder));

    for ix=1:sum(sum(TTmode2==1)),
        eEz{iorder,2}(ix)=1/(SEz{iorder+1}(ix)-SEz{iorder}(ix));
        eHx{iorder,2}(ix)=1/(SHx{iorder+1}(ix)-SHx{iorder}(ix));
        eHy{iorder,2}(ix)=1/(SHy{iorder+1}(ix)-SHy{iorder}(ix));
    end;

end;
for ip=3:N_order,
    for iorder=1:N_order-ip+1,
        e{ip}(iorder)=e{ip-2}(iorder+1)+1/(e{ip-1}(iorder+1)-e{ip-1}(iorder));    

        for ix=1:sum(sum(TTmode2==1)),
            eEz{iorder,ip}(ix)=eEz{iorder+1,ip-2}(ix)+1/(eEz{iorder+1,ip-1}(ix)-eEz{iorder,ip-1}(ix));    
            eHx{iorder,ip}(ix)=eHx{iorder+1,ip-2}(ix)+1/(eHx{iorder+1,ip-1}(ix)-eHx{iorder,ip-1}(ix));    
            eHy{iorder,ip}(ix)=eHy{iorder+1,ip-2}(ix)+1/(eHy{iorder+1,ip-1}(ix)-eHy{iorder,ip-1}(ix));    
        end;

    end;
end;
%prediction(ialpha)=log10((abs(e{N_order})).^2)+2.369;%log10((abs(e{N_order})).^2)+2.369;

%prediction(ialpha)=1/(16*Pmode*Pinc)*(0.05*abs(sum((eEz{1,N_order}).'.*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((eHx{1,N_order}).')))).^2;
prediction(ialpha)=1/(16*Pmode*Pinc)*(0.05*abs(sum((eEz{N_order,1}).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((eHx{N_order,1}))))).^2;

%prediction(ialpha)=1/(16*Pmode*Pinc)*(0.05*abs(sum((SEz{N_order}).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((SHx{N_order}))))).^2;
%prediction(ialpha)=Pmode*Pinc*(abs(e{N_order})).^2;


%prediction(ialpha)=1/(4*sqrt(Pmode*Pinc))*(0.05*(sum((SEz{N_order}).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((SHx{N_order})))));
%prediction(ialpha)=1/(4*sqrt(Pmode*Pinc))*(0.05*(sum((eEz{1,N_order}).'.*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((eHx{1,N_order}).'))));

% eps_str_brute=eps_str_old;
% eps_str_brute(TTopt==1)=eps_str_brute(TTopt==1)+dFOM(TTopt==1).*(alphav(ialpha)).*facmin;%0.005
% Mz = zeros(N);
% Mz(TTmode1==1)=(modeTE0re+i.*modeTE0im);
% [Ez_for_brute, Hx_for_brute, Hy_for_brute] = solveTM(L0, wvlen, xrange, yrange, eps_str_brute, Mz, Npml);
% obj_brute(ialpha)=1/(16*Pmode*Pinc)*(0.05*abs(sum((Ez_for_brute(TTmode2==1)).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((Hx_for_brute(TTmode2==1)))))).^2;

end;

%%en brut
% UUUint=0:0.001:1;
% VVVint=interp1(alphav,full(obj_brute),UUUint,'spline');
% WWW(iiter)=UUUint(find(VVVint==max(VVVint)));

%if iiter==4,return,end;

% %%%en prediction
UUUint=1e-6:1e-6:1;
VVVint=interp1(alphav,full(prediction),UUUint,'spline');
WW=find(VVVint==max(VVVint));
WWW(iiter)=UUUint(WW(1));


%prediction(97)
if iiter==1,
    alphamax(iiter)=WWW(iiter);%alphafix;%alphafix;%alphafix;%0.2;%alphafix;%WWW(iiter);%0.5;% 1.8239e-10 0.0118 0.0459 0.0981 0.1620 0.2298 0.2934 0.3456 0.3813 0.3974 0.3935
%0.6130 + 0.1501i 0.6130 + 0.1501i 0.6130 + 0.1501i 0.6130 + 0.1501i 0.6130 + 0.1501i
else,
    alphamax(iiter)=WWW(iiter);%alphafix;%alphafix;%alphafix;%0.2;%alphafix;%WWW(iiter);
end%alphav(find(prediction==max(prediction))),end;
%prediction=log10((abs(e{N_order})).^2)+2.369
figure(75),plot(alphav,prediction)

if 0,%iiter==2,
    return,
end;

for iorder=2:N_order,
SEz{iorder}=SEz{iorder-1}+(alphamax(iiter)).^(iorder-1)*dEz_for{iorder}(TTmode2==1);%dEz_for{1}(TTmode2==1);
SHx{iorder}=SHx{iorder-1}+(alphamax(iiter)).^(iorder-1)*dHx_for{iorder}(TTmode2==1);
SHy{iorder}=SHy{iorder-1}+(alphamax(iiter)).^(iorder-1)*dHy_for{iorder}(TTmode2==1);
end;
for iorder=1:N_order,
SEz{iorder}=full(SEz{iorder});
SHx{iorder}=full(SHx{iorder});
SHy{iorder}=full(SHy{iorder});
end;
for iorder=1:N_order,
eEz{iorder,1}=SEz{iorder};
eHx{iorder,1}=SHx{iorder};
eHy{iorder,1}=SHy{iorder};
end;
for iorder=1:N_order-1,
for ix=1:sum(sum(TTmode2==1)),
eEz{iorder,2}(ix)=1/(SEz{iorder+1}(ix)-SEz{iorder}(ix));
eHx{iorder,2}(ix)=1/(SHx{iorder+1}(ix)-SHx{iorder}(ix));
eHy{iorder,2}(ix)=1/(SHy{iorder+1}(ix)-SHy{iorder}(ix));
end;
end;
for ip=3:N_order,
    for iorder=1:N_order-ip+1,
for ix=1:sum(sum(TTmode2==1)),
eEz{iorder,ip}(ix)=eEz{iorder+1,ip-2}(ix)+1/(eEz{iorder+1,ip-1}(ix)-eEz{iorder,ip-1}(ix));    
eHx{iorder,ip}(ix)=eHx{iorder+1,ip-2}(ix)+1/(eHx{iorder+1,ip-1}(ix)-eHx{iorder,ip-1}(ix));    
eHy{iorder,ip}(ix)=eHy{iorder+1,ip-2}(ix)+1/(eHy{iorder+1,ip-1}(ix)-eHy{iorder,ip-1}(ix));    
end;
    end;
end;


for t=2:N_order,
eEz{1,t}=(eEz{1,t}).';
eHx{1,t}=(eHx{1,t}).';
eHy{1,t}=(eHy{1,t}).';
end;

% figure(70),plot(XX(TTmode2==1),real(eEz{1,N_order}),'r'),hold on,for t=1:N_order,plot(XX(TTmode2==1),real(SEz{t}),'k'),end;hold off,xlabel('X(µm)')
% for t=1:N_order,
% Gnu(t)=1/(4*sqrt(Pmode*Pinc))*(0.05*(sum((SEz{t}).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((SHx{t})))));
% Gnu2(t)=1/(4*sqrt(Pmode*Pinc))*(0.05*(sum((eEz{1,t}).*conj(Hx_mode(TTmode2==1))+conj(Ez_mode(TTmode2==1)).*((eHx{1,t})))));
% end;
% Gnu2(1)=Gnu(1);
% figure(71),hold,plot(0:N_order-1,real(Gnu),'k'),plot(0:2:N_order-1,real(Gnu2(1:2:N_order)),'r')
% figure(72),hold,plot(0:N_order-1,imag(Gnu),'k'),plot(0:2:N_order-1,imag(Gnu2(1:2:N_order)),'r')
% figure(73),hold,plot(0:N_order-1,log10(abs(Gnu-(0.61302 + 0.15011i))./abs(0.61302 + 0.15011i)),'k'),plot(0:2:N_order-1,log10(abs(Gnu2(1:2:N_order)-(0.61302 + 0.15011i))./abs(0.61302 + 0.15011i)),'r')
% figure(74),hold,plot(0:N_order-1,log10(abs(Gnu-(Gnu2(end)))./abs(Gnu2(end))),'k'),plot(0:2:N_order-1,log10(abs(Gnu2(1:2:N_order)-(Gnu2(end)))./abs(Gnu2(end))),'r')
% 

%eps_str(TTopt==1)=eps_str_old(TTopt==1)+dFOM(TTopt==1).*(1./abs(dEz_for{iorder}(II1b,II2b))).*0.1*facmin;%0.005
eps_str(TTopt==1)=eps_str(TTopt==1)+dFOM(TTopt==1).*(alphamax(iiter)).*facmin;%0.005




% MMC=[reshape(dEz_for{1},141*141,1)-reshape(dEz_for{2},141*141,1)];
% MMC2=[reshape(dEz_for{1},141*141,1)];
% for iorder=2:N_order-1,
% MMC=[MMC reshape(dEz_for{iorder},141*141,1)-reshape(dEz_for{iorder+1},141*141,1)];
% MMC2=[MMC2 reshape(dEz_for{iorder},141*141,1)];
% end;
% FMC=[reshape(dEz_for{1},141*141,1)];
% alpha=MMC\FMC;
% Ez_proj2{iiter}=reshape(MMC2*alpha,141,141);

%predictedObjective=log10((abs(Ez_proj{iiter}(II1b,II2b)))^2)+2.369




end;

if newbinaire,%iiter==70,
    
    VV=eps_str;
    VV(eps_str>(ff*1.5^2+(1-ff)*2.5^2))=1;
    eps_str(eps_str>(ff*1.5^2+(1-ff)*2.5^2))=2.5^2;
    
    VV(eps_str<=(ff*1.5^2+(1-ff)*2.5^2))=-1;
    eps_str(eps_str<=(ff*1.5^2+(1-ff)*2.5^2))=1.5^2;
    
    VV2=VV;
    figure(10),pcolor(VV.'),shading('interp'),colorbar
    
    for ix=1:size(VV2,1),
        for iy=1:size(VV2,2),
            
            trouve=0;
            t=1;
            
      while trouve==0
          if any(any([sign(VV(max(1,ix-t):min(size(VV,1),ix+t),max(1,iy-t):min(size(VV,2),iy+t)))==-sign(VV(ix,iy))])),
              VV2(ix,iy)=sign(VV(ix,iy)).*t;
              trouve=1;
          end;
          t=t+1;
          end;
              
            %VV2(ix,iy)=VV(ix,iy)+VV(ix,iy+1)+VV(ix,iy-1)+VV(ix+1,iy)+VV(ix+1,iy+1)+VV(ix+1,iy-1)+VV(ix-1,iy)+VV(ix-1,iy+1)+VV(ix-1,iy-1);
      end;
  
  end;
    figure(11),pcolor(VV2.'),shading('interp'),colorbar
    
   newbinaire=0;
   disp('on a construit le ls')
        end;
   
  

if binaire,%iiter>70,

    
% if max(dFOM(TTopt==1).*(2.5^2-eps_str(TTopt==1)))>-min(dFOM(TTopt==1).*(eps_str(TTopt==1)-1)),signdf=1;else,signdf=0;end;
% if signdf,
% eps_str(and(TTopt==1,(dFOM.*(2.5^2-eps_str))==max(dFOM(TTopt==1).*(2.5^2-eps_str(TTopt==1)))*1))=2.5^2;
% else,
% eps_str(and(TTopt==1,(dFOM.*(eps_str-1.5^2))==min(dFOM(TTopt==1).*(eps_str(TTopt==1)-1.5^2))*1))=1.5^2;
% end;

%if max(dFOMp(TTopt==1).*(2.5^2-eps_str(TTopt==1)))>max(dFOMn(TTopt==1).*(eps_str(TTopt==1)-1.5^2)),signdf=1;else,signdf=0;end;
%if signdf,
%eps_str(and(TTopt==1,(dFOMp.*(2.5^2-eps_str))==max(dFOMp(TTopt==1).*(2.5^2-eps_str(TTopt==1)))*1))=2.5^2;
%else,
%eps_str(and(TTopt==1,(dFOMn.*(eps_str-1.5^2))==min(dFOMn(TTopt==1).*(eps_str(TTopt==1)-1.5^2))*1))=1.5^2;
%end;


%chang=0;
%while chang==0,
%VV2b=VV2+0.00001.*dFOM;
%if any(any(VV2b.*VV2<=0)),chang=1;end;
%VV2=VV2b;
%end;

%VV2=VV2+0.05.*dFOM;


%levelset
dt=0.2*1+0.5*0;%0.5./(1+cptbinaire/50);%0.1

normalisls=max(abs(dFOM(TTopt==1)));

VV2b=VV2;
for ix=2:size(VV2,1)-1,
        for iy=2:size(VV2,2)-1,
    
            if dFOM(ix,iy)<0,            
            VV2(ix,iy)=VV2b(ix,iy)+dFOM(ix,iy)./normalisls.*dt.*sqrt((min(VV2b(ix+1,iy)-VV2b(ix,iy),0))^2+(max(VV2b(ix,iy)-VV2b(ix-1,iy),0))^2+(min(VV2b(ix,iy+1)-VV2b(ix,iy),0))^2+(max(VV2b(ix,iy)-VV2b(ix,iy-1),0))^2);
            else,
            VV2(ix,iy)=VV2b(ix,iy)+dFOM(ix,iy)./normalisls.*dt.*sqrt((max(VV2b(ix+1,iy)-VV2b(ix,iy),0))^2+(min(VV2b(ix,iy)-VV2b(ix-1,iy),0))^2+(max(VV2b(ix,iy+1)-VV2b(ix,iy),0))^2+(min(VV2b(ix,iy)-VV2b(ix,iy-1),0))^2);
                
            end;
        end;
end;
eps_str((TTopt.*VV2)>0)=2.5^2;
eps_str((TTopt.*VV2)<0)=1.5^2;

%cptbinaire=cptbinaire+1;
% 
% if cptbinaire>(1000),
% irandomfini=1;
%  clear Objective
% end;

end;

%  XXX=XX(TTopt==1);
%  YYY=YY(TTopt==1);
%  dFFOM=dFOM(TTopt==1);
%  
%  alpha=0.01./max(abs(dFOM(TTopt==1)));
%  
%  for iX=1:numel(XXX),
%      
%      if dFFOM(iX)>0,
%      eps_str(and(TTopt==1,sqrt((XX-XXX(iX)).^2+(YY-YYY(iX)).^2)<alpha.*(dFFOM(iX)-max(abs(dFOM(TTopt==1)))*0.9)))=4;
%      else,
%          
%      eps_str(and(TTopt==1,sqrt((XX-XXX(iX)).^2+(YY-YYY(iX)).^2)<alpha.*(-dFFOM(iX)-max(abs(dFOM(TTopt==1)))*0.9)))=1;
%      end;
%      
%  end;

%figure
%pcolor(XX,YY,eps_str_old);shading('interp'),title('epsilon old'),caxis([2.25 6.25]),xlabel('X(um)'),ylabel('Z(um)'),colorbar,axis([min(min(XX0b)) max(max(XX0b)) min(min(YY0b)) max(max(YY0b))]);


%figure
%pcolor(XX,YY,eps_str);shading('interp'),title('epsilon'),caxis([2.25 6.25]),xlabel('X(um)'),ylabel('Z(um)'),colorbar,axis([min(min(XX0b)) max(max(XX0b)) min(min(YY0b)) max(max(YY0b))]);
frame = getframe(figure(4));
im4{iiter} = frame2im(frame);

drawnow


iiter=iiter+1;


end;



end;




%% Visualize the solution.
%visabs(Ez, xrange, yrange);title('Ez')
%pcolor(real(Ez));shading('interp'),colorbar

%%
%figure;
%visabs(Hx, xrange, yrange - 0.5*dL(2));title('Ex')
%pcolor(abs(Ex));shading('interp'),colorbar

%figure;
%visabs(Hy, xrange - 0.5*dL(1), yrange);title('Ey')
%pcolor(abs(Ey));shading('interp'),colorbar

%figure,pcolor(linspace(xrange(1),xrange(2),N(1)),linspace(yrange(1),yrange(2),N(2)),eps_air),shading('interp'),colorbar

%% Show the movie of the oscillating field.
% figure;
% moviereal(Hz, xrange, yrange)


figure();
plot(alphamax,'LineWidth', 2);
xlabel('iteration')
ylabel('alpha')
xlim([0,400])
