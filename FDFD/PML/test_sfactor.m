%% sfactor test cases
%% CREATE SFACTOR ALREADY ANTICIPATES THAT YOU HAVE CREATED A CONSISTENT DOMAIN PML
%%===============0 length PML========================00
L0 = 1e-6
wrange = [-1,1]
s = 'f'
eps0 = 8.85*10^-12*L0;
mu0 = 4*pi*10^-7*L0;
c = 1/sqrt(mu0*eps0); c0 = c;

Nw = 200;
Nw_pml = 20;
Nw = Nw; %this is always true in the actual code;
M = Nw*Nw;
L0 = 10^-6; wvlen = 1;
omega = 2*pi*c0/(wvlen);
lnR = -16;

sfactor_array = create_sfactor(wrange, s, omega, eps0, mu0, Nw, Nw_pml);
figure;
plot(abs(sfactor_array));

%%===================Excessive Length PML======================

sfactor_array_mine = create_sfactor_mine(wrange, s, omega, eps0, mu0, Nw, Nw_pml);
figure;
plot(abs(sfactor_array_mine));

%%================Compare to hw solutoin==================%
figure;
plot(abs(sfactor_array));
hold on;
plot(abs(sfactor_array_mine))

%% create a bunch of sfactors
sfactor_f= create_sfactor(wrange, 'f', omega, eps0, mu0, 200, Nw_pml);
sfactor_b= create_sfactor(wrange, 'b', omega, eps0, mu0, 200, Nw_pml);
figure();


plot(imag(sfactor_f),'.');
hold on;
plot(imag(sfactor_b),'.');


%% actual code
lnR = -12;  % R: target reflection coefficient for normal incidence

m = 3.5;% degree of polynomial grading
%% Output Parameter
% sfactor_array: 1D array with Nw elements containing PML s-factors for Dws

eta0 = sqrt(mu0/eps0);  % vacuum impedance

w_array = linspace(wrange(1), wrange(2), Nw+1);

loc_pml = [w_array(1 + Nw_pml), w_array(end - Nw_pml)]; % specifies where the pml begins on each side
d_pml = abs(wrange - loc_pml); % pml thickness

%% what happens when we have a 0 in the denominator
sigma_max = -(m+1)*lnR/(2*eta0) ./ d_pml; %usually the pml is the same thickness on both sides

%% forward or backward...idk what this is requiring
if s == 'b' % ws is 1 unit smaller than w_array
    ws = w_array(1:end-1);
else  % s == 'f'
    assert(s == 'f');
    ws = (w_array(1:end-1) + w_array(2:end)) / 2;
end

ind_pml = {ws < loc_pml(1), ws > loc_pml(2)};  % {} means cell array

sfactor_array = ones(1, Nw);
for n = 1:2
    sfactor = @(L) 1 - 1i * sigma_max(n)/(omega*eps0) * (L/d_pml(n)).^m;
    sfactor_array(ind_pml{n}) = sfactor(abs(loc_pml(n) - ws(ind_pml{n})));
end

%% generate the array
%% these are row vectors
sfactor_f= create_sfactor(wrange, 'f', omega, eps0, mu0, 200, Nw_pml);
sfactor_b= create_sfactor(wrange, 'b', omega, eps0, mu0, 200, Nw_pml);

[Sxf, Syf] = ndgrid(sfactor_f,sfactor_f);
[Sxb, Syb] = ndgrid(sfactor_b,sfactor_b);

Sxf=spdiags(Sxf(:),0,M,M);
Sxb=spdiags(Sxb(:),0,M,M);
Syf=spdiags(Syf(:),0,M,M);
Syb=spdiags(Syb(:),0,M,M);

%% what should the structure of this array really be?
validate = imag(diag(Sxf));
figure(); plot(validate(1:201))

validate = imag(diag(Syf));
figure(); plot(validate(1:301))
