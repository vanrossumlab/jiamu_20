% using matlab to solve PDE for voltage on a sphere given current injection
% variant w cos(th)-transform
close all
clear all

global tau rho d rm ri rholam2 z_a Iext t_on t_off
% Unit of parameters:
% Length: Micrometer 
% Resistor: Ohm
% Capacitor: Farad

ri  = 1;    % internal resistivity
rm  = 1;    % leak per area
d   = 1;    % shell thickness
Iext= 1;    % total stimulus current, applied between theta=0..theta_a
th_a= 0.1;

t_on=0;
t_off=0.02;
Iext=Iext / (t_off-t_on);

lam= sqrt(rm*d/ri);

tau=1; % we assume that cm is set so that tau=rm*cm =1 

nz      = 500;   % # of points used for z
z_a     = cos(th_a);
z_min   = cos(pi);
z_ar    = linspace(z_min, z_a, nz);

nt      = 3000;    % # of points used for time
tmin    = 0;
tmax    = 5;
t_ar    = linspace(tmin,tmax,nt);

m=0; % no additional symmetry

rho_ar = [0.1 0.5 1 2]
lam= sqrt(rm*d/ri);

res=t_ar';
for rho = rho_ar
    rholam2= (rho/lam)
    sol = pdepe(m,@pdefun,@icfun,@bcfun,z_ar,t_ar);
    vpi_vs_t = sol(:,1);

    figure(1)
    hold on
    plot(t_ar,vpi_vs_t)
    drawnow

    figure(2) %nromalized
    hold on
    plot(t_ar,vpi_vs_t/max(vpi_vs_t))
    drawnow
    res=[res vpi_vs_t];
end
%save sphere_pulse.dat res -ascii


function [c,f,s] = pdefun(z,t,u,dudz)
    global tau rholam2
    
    c = tau*rholam2;
    
    f = (1-z^2)*dudz;
    
    s = -rholam2*u;
end

function [pL,qL,pR,qR] = bcfun(zL,uL,zR,uR,t)
    global Iext ri d z_a t_on t_off
    
    % BCs defined as : p+q*f(z,t,u,dudz)=0
    % here f=(1-z^2)*dudz=sin^2*dudz
    % V'(z_a) = -Iext *ri/(2*pi*d)/sin^2(za)

    pL = 0; % note Left is th=pi,z=-1
    
    qL = 1; % should be automatic
    
    pR = -Iext*ri/(2*pi*d)*(t<t_off).*(t>t_on);
    
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
    
    % dz=abs(z_ar(2)-z_ar(1));
    % ileak=sum(u)*dz; 
    ileak = trapz(z_ar,u); % or integral-command
    ileak= ileak*2*pi*rho^2/rm;
end

