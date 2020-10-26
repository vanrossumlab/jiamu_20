close all
clear all
% model parameters
global tau rho d rm ri rholam2 z_a z_stim
global vclamp
global stim_amp stim_z_min stim_z_max stim_t_on stim_t_off th_a


ri  = 1;    % internal resistivity
rm  = 1;    % leak per area
d   = 1;    % shell thickness
rho = 1;    % diameter of sphere
th_a= 0.1;

lam = sqrt(rm*d/ri);
plotlist=[3];
tau=1; % we assume that cm is set so that tau=rm*cm =1 

nz      = 2000;   % # of points used for z
z_a     = cos(th_a); % size  of clamp electrode
th_stim = th_a;
z_stim   = cos(pi-th_stim); % size of current input at bottom
z_ar    = linspace(z_stim, z_a, nz);
dz      = (z_a-z_stim)/nz;

nt      = 1000;    % # of points used for time
tmin    = 0;
tmax    = 10;
t_ar    = linspace(tmin,tmax,nt);

dt=tmax/nt;
m=0;

% stimulus parameters
vclamp = 0;
stim_amp = 1; % current density
stim_t_on = 5;
stim_t_off = tmax;
it_stim= stim_t_on/dt;
t_ar_stim=t_ar(it_stim:end)-stim_t_on;

res=t_ar_stim';
res2=[];

%for  rho =[0.1 0.5 1 2 10]% for sample run
for  rho=0.1:0.1:10 % for stats run
    rho_over_lam = (rho/lam) % this variable determines how compact the sphere is 
    rholam2=rho_over_lam^2;
    
    %solve partial differential equation
    sol = pdepe(m,@pdefun,@icfun,@bcfun,z_ar,t_ar);

     if (find(plotlist==1))
        figure(1)
        surf(sol);
        %u = sol(:,:,1);
        %surf(z_ar,t_ar,u);
        xlabel('z')
        ylabel('t')
        zlabel('V')
     end   
     dv=-(sol(:,end-1)-sol(:,end))/dz; % -V'(x=0,t)
     d_th=pi/nz;
     iclamp = dv*(2*pi*d)*(sin(th_a))^2/ri;
    
    ileak=zeros(1,nt);
    istim=zeros(1,nt);
    for it=1:nt
        ileak(it) = -trapz(z_ar,sol(it,:)); % need minus..?
        istim(it) = stim_fun(it*dt);
       % istim(it)=0;
       % for z=z_ar
            %istim(it) = istim(it)+stim_fun(it*dt)*dz;
       % end
    end    
    ileak=ileak*2*pi*rho^2/rm;
   
    % check current conservation, 
    % for now without transient capacitive current
     
    if (find(plotlist==2))
        figure(2);
        plot(t_ar(it_trans:end),iclamp)
        xlabel("time"); ylabel("current");
        hold on
        plot(t_ar(it_trans:end),ileak)
        plot(t_ar(it_trans:end),istim)
        plot(t_ar(it_trans:end),isum)

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
        figure(3)
        plot(t_ar_stim,neg_iclamp);
        drawnow
        hold on
    end
    res=[res neg_iclamp];
    res2=[res2 [rho t50 iclamp_ss]'];
end
%save sphere_vclamp.dat res -ascii
res2=res2';
save sphere_vclamp_stats.dat res2 -ascii


function stim=stim_fun(t)
% returns current density
    global stim_amp stim_t_on stim_t_off
    stim=stim_amp*(t>stim_t_on)*(t< stim_t_off);
end    

function [c,f,s] = pdefun(z,t,u,dudz)
    global tau rholam2 rm rho z_a th_a
    
    c = tau*rholam2;
    
    f = (1-z^2)*dudz;
    
    s = -rholam2*u;
    
   % s = s + rholam2*stim_fun(t,z)*rm/(2*pi*rho^2);
end

function [pL,qL,pR,qR] = bcfun(zL,uL,zR,uR,t)
    global vclamp d ri th_a th_stim
  % BCs defined as : p+q*f(z,t,u,dudz)=0, where f is flux

    % note Left is th=pi,z=-1  , Here p=-f;
    pL = stim_fun(t)/(2*pi*d)/ri;
    qL = 1 ;
    
    % northpole is clamped
    pR = uR-vclamp;
    qR = 0; 
end


function u0 = icfun(z)  % initial condition
    u0 = 0;
end

