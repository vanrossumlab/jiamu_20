% using matlab to solve PDE for voltage on a sphere given current injection
% variant w cos(th)-transform
close all
clear 
% clc

global tau rho d rm ri rholam2 z_a Iext


ri  = 100;    % internal resistivity
rm  = 10^8;    % leak per area
d   = 0.3;    % shell thickness
rho = 600;    % diameter of sphere
Iext= 10;    % total stimulus current, applied between theta=0..theta_a
th_a= 0.1;

lam= sqrt(rm*d/ri);

rholam2=1;% (rho/lam)^2 % this variable determines how compact the sphere is 

tau=1; % we assume that cm is set so that tau=rm*cm =1 

nz      = 500;   % # of points used for z
z_a     = cos(th_a);
z_min   = cos(pi);
z_ar    = linspace(z_min, z_a, nz);

nt      = 500;    % # of points used for time
tmin    = 0;
tmax    = 7;
t_ar    = linspace(tmin,tmax,nt);

m=0; % no additional symmetry
sol = pdepe(m,@pdefun,@icfun,@bcfun,z_ar,t_ar);
% pdepe: Solve initial-boundary value problems for parabolic-elliptic PDEs in 1-D

% some sanity checks
uss= sol(nt,:);
Vpeak_ss = uss(1)

Asphere=4*pi*rho^2;
uss_sc = Iext/Asphere*rm % theory:steady state single comp. small th_a
totalLeak = total_leak(z_ar,uss);
fprintf("error in steady state(1-Iext/Ileak) = %g\n", 1-Iext/totalLeak);

% make charging time plot
    t_ar_plot=t_ar(t_ar<100);
    nplot=length(t_ar_plot);
    vpi_vs_t = sol(:,1);
    va_vs_t = sol(:,nz);

    vpi_vs_t =vpi_vs_t/max(va_vs_t);
    va_vs_t =va_vs_t /max(va_vs_t);

    figure
    
    hold on
    plot(t_ar_plot,va_vs_t(1:nplot),'DisplayName','Sphere')

    exp_model = (1-exp(-t_ar_plot/tau));
    plot(t_ar_plot, exp_model,'DisplayName','Cable')
    erf_model = erf(sqrt(t_ar_plot/tau)); % inf cable
    plot(t_ar_plot, erf_model,'DisplayName','Single compartment')
    legend
    
%     csvwrite('large_sphere_charge.csv',va_vs_t);
   
function [c,f,s] = pdefun(z,t,u,dudz)
    global tau rholam2
    
    c = tau*rholam2;
    
    f = (1-z^2)*dudz;
    
    s = -rholam2*u;
end

function [pL,qL,pR,qR] = bcfun(zL,uL,zR,uR,t)
    global Iext ri d z_a
  
    % note Left is th=pi,z=-1
    
    % BCs defined as : p+q*f(z,t,u,dudz)=0
    % here f=(1-z^2)*dudz=sin^2*dudz
    % V'(z_a) = -Iext *ri/(2*pi*d)/sin^2(za)

    pL = 0;
    
    qL = 1; % should be automatic
    
    pR = -Iext*ri/(2*pi*d);
    
    qR = 1; 
end


function u0 = icfun(z) 
    % initial condition
    u0 = 0;
end

function ileak=total_leak(z_ar,u)
    global rho rm
    % given solution u(z), calculate total leak
    % the result should in steady state equal total injected current
    % assumes equal z spacing
     
    % for single comp this should equal V/R=uss(1)/rm*Asphere
    
    dz=abs(z_ar(2)-z_ar(1));
    
    ileak=sum(u)*dz; 
    
    ileak = trapz(z_ar,u); % or integral-command
    
    ileak=ileak*2*pi*rho^2/rm;
end

