% examines voltage clamp in finite cable.
% specifically this script writes cable_vclamp.dat that contains
% clamp current for a number of cable lengths.
% The input is always at the end of the cable, and implemented through a
% boundary condition.

close all
clear all
%declare model parameters
global C tau rm ra d vclamp
global stim_amp stim_t_on stim_t_off

C=1;
d=4; % makes lambda =1

ra=1;
rm=1;
tau=1;
lam2_mvr=(d*rm/ra/4)

tmax=10;

vclamp= 0
nt=1000;
t_ar = linspace(0,tmax,nt);
dt=tmax/nt;

plotlist=[3];

% stimulus parameters
stim_amp = 1; % current density
stim_t_on = 5;
stim_t_off = tmax;
it_stim= stim_t_on/dt;
t_ar_stim=t_ar(it_stim:end)-stim_t_on;

res=t_ar_stim';
res2=[];

%for  L=[0.1 0.5 1 2 10] % for samples
for  L=2.7*[0.1 0.5 1 2 10] % for samples
%for  L=0.1:0.1:10 % for stats run
    L
    lam2_mvr=(d*rm/ra/4);
    
    L_over_lam=L/sqrt(lam2_mvr)
    
    m=0;
    nx=500;
    x_ar = linspace(0.0,L,nx);
    dx= L/nx;
    
    %solve partial differential equation
    sol = pdepe(m, @linear_cable_pde_eqn, @linear_cable_pde_initial, @linear_cable_pde_bc, x_ar,t_ar);
    
    if (find(plotlist==1))
        figure
        surf(sol);
        
        %u = sol(:,:,1);
        %surf(x_ar,t_ar,u);
        xlabel('x')
        ylabel('t')
        zlabel('V')
    end
    
    dv=-(sol(:,2)-sol(:,1))/dx; % -V'(x=0)
    iclamp = dv*pi*d^2/(4*ra);
    
    for it=1:nt
        ileak(it) = trapz(x_ar,sol(it,:));
        istim(it)=0;
        for x=x_ar
            istim(it) = istim(it)+stim_fun(it*dt)*dx;
        end
    end
    ileak = ileak*pi*d/rm;
    
    isum = iclamp+ileak+istim; % should sum to zero in SS ( Icap not calculated)
    
    
    if (find(plotlist==2))
        figure(2);
        plot(iclamp)
        xlabel("time"); ylabel("current");
        hold on
        % check current conservation, for now without transient capacitive current
        plot(ileak)
        plot(istim)
        plot(isum)
        legend("clamp","leak","stim","sum");
    end
    
    % note for vclamp != 0, we should subtract iclampSS
    iclamp = iclamp(it_stim:end);
    neg_iclamp= -iclamp;
    ileak = -ileak(it_stim:end)';
    istim = istim(it_stim:end)';
    
    iclamp_ss =max(neg_iclamp)
    t50 =dt*find(neg_iclamp>iclamp_ss/2,1)
    
    if (find(plotlist==3))
        plot(t_ar_stim,neg_iclamp);
        drawnow
        hold on
    end
    res=[res neg_iclamp];
    res2=[res2 [L t50 iclamp_ss]'];
end
save cable_vclamp_scaled.dat res -ascii
res2=res2';
%save cable_vclamp_stats.dat res2 -ascii


function [c,b,s] = linear_cable_pde_eqn(x,t,u,DuDx)
global rm ra d tau vclamp;
% c* dudt = x^-m d(x^m f)/dx +s
% here m=0; s,f,c functions of (x,t,u,dudx)
% lam=(rm/ra)^0.5; like sphere
lam2_mvr=(d*rm/ra/4); % standard definition.

c = tau;
b = lam2_mvr*DuDx;
s = -u;

%input for vclamp. % Now moved to BCfun
% s = s + stim_fun(t,x)*rm/pi/d;
end

function stim=stim_fun(t)
% returns current density, takes a lot of CPU time!
global stim_amp stim_t_on stim_t_off
stim=stim_amp*(t>stim_t_on)*(t< stim_t_off);
end

%set initial condition
function value = linear_cable_pde_initial(x)
value = 0;
end

%set boundary conditions
function [pl,ql,pr,qr] = linear_cable_pde_bc(zL,uL,zR,uR,t)
global ra d vclamp rm
% BCs defined as : p+q*f(z,t,u,dudz)=0, where f is flux

pl = uL-vclamp;
ql = 0;
pr = -stim_fun(t)/(pi*d*rm);
qr = 1;
end
