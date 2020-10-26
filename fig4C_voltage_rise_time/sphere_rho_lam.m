% using matlab to solve PDE for voltage on a sphere given current injection
% variant w cos(th)-transform
close all
clear all

% Unit of parameters:
% Length: Micrometer 
% Resistor: Ohm
% Capacitor: Farad

ri  = 1;    % internal resistivity
rm  = 1;    % leak per area
d   = 1;    % shell thickness
Iext= 1;    % total stimulus current, applied between theta=0..theta_a

lam= sqrt(rm*d/ri)

tau=100; % we assume that cm is set so that tau=rm*cm =1 

nz      = 2000;  % # of points used for z

nt      = 16000;    % # of points used for time
tmax    = 10*tau;
t_ar    = linspace(0,tmax,nt);
dt= tmax/nt;


m=0; % no additional symmetry
t_ar_plot=t_ar;
% get the maxima and minima of t50
exp_model = (1-exp(-t_ar_plot/tau));
t50exp = t_ar(find(exp_model>0.5*max(exp_model),1)); % 50% time;
erf_model = erf(sqrt(t_ar_plot/tau));
t50erf = t_ar(find(erf_model>0.5*max(erf_model),1)); % 50% time;


simcase=2
switch(simcase)
case(1) % t50 across number of radii and electrode angle

    rho_ar = [lam/10, lam/5, lam/2, lam, 2*lam 10*lam 20*lam];
    tha_ar = 0.005:0.01:0.9; %t came in near perfect condition, complete with power & data cable + toner. The price was fantastic and the printer has been working like a trojan since it was installed. We would definitely buy this printer again, great speed and quality.

    nrho=length(rho_ar);
    nth=length(tha_ar);


    t50=zeros(nrho,nth);
    vss=zeros(nrho,nth);
    %for j=1:nrho
    parfor j=1:nrho
        rho = rho_ar(j);
        rholam2= (rho/lam)^2        
        
        pprr=-Iext*ri/(2*pi*d);
        for i=1:nth
            th_a = tha_ar(i);
            z_a     = cos(th_a);
            z_min   = cos(pi);
            z_ar    = linspace(z_min, z_a, nz);
            sol = pdepe(m,@(z,t,u,dudz)pdefun(z,t,u,dudz, tau,rholam2),@icfun,@(zL,uL,zR,uR,t)bcfun(zL,uL,zR,uR,t,pprr),z_ar,t_ar);
            va_vs_t = sol(:,nz); % v(t,th=theta_a)
            % for debug, only works with regular for, not parfor
            %figure(3); plot(t_ar,va_vs_t); hold oncloase all
            toff=10*tau;
            vss=sol(ceil(toff/dt-1),:);
            
            sol_plot=downsample(downsample(sol,ceil(nt/100))',ceil(nz/100))';
            figure(1)
                surf(sol_plot)
            
            figure(2)
                plot(t_ar,va_vs_t)
                ylabel('V at pipette')
                xlabel('time')
                
            figure(3)
                plot(vss);
                xlabel('position  as cos(theta)')
                ylabel('voltage at peak')
            
            va_vs_t = va_vs_t /max(va_vs_t);
            it=find(va_vs_t>0.5*max(va_vs_t),1);
            if (it<10)
                printf("reaching time res, increase nt\n")
            end        
            t50(j,i) = t_ar(it); % 50% time
            vss(j,i) = max(sol(:,nz)); % input resistance
        end
    end
    figure(1)
        set(gcf,'position',[0,0,400,300])
        t50_exp= t50exp * ones(length(tha_ar),1);
        t50_erf= t50erf * ones(length(tha_ar),1);
        plot(tha_ar,t50,'-','LineWidth',1)
        hold on
        plot(tha_ar,t50_erf,'LineWidth',1);
        plot(tha_ar,t50_exp,'LineWidth',1);

        legend(string(rho_ar))
        xlabel('\theta_a')
        ylabel('Half maximum time (\tau)')
        ylim([0 0.75])
        annotation('textbox', [0.45, 0.4, 0.18, 0.05], 'String', 'Infinite Cable','EdgeColor','white')
        annotation('textbox', [0.45, 0.8, 0.4, 0.05], 'String', 'Single Compartment','EdgeColor','white')
        savefig('sphere_t50_rho_lam.fig')
        saveas(gcf,'sphere_t50_rho_lam.png')

    figure(2)    
        plot(tha_ar,vss,'-','LineWidth',1)
        legend(string(rho_ar))
        xlabel('\theta_a')
        ylabel('Vss')

    t50=t50';
    t50p=[tha_ar' t50 t50_exp t50_erf];

    csvwrite('t50_sphere_rho_lam.csv',t50p)

    csvwrite('vss_sphere_rho_lam.csv',vss)

case(2)
% constant electrode SIZE
    rho_ar = lam*0.01*[1.1:2:20];

    nrho=length(rho_ar);
    
    s_a_ar=[1e-3 0.005, 0.01];
    nth=length(s_a_ar);
    t50=zeros(nrho,nth);
    vss=zeros(nrho,nth);
    for i=1:nth
        s_a=s_a_ar(i);
        for j=1:nrho
            rho = rho_ar(j);

            if (s_a>= rho)
               fprintf(1,'s_a = %g,   rho= %g \n', s_a,rho)
               %error('stop s_a>= rho')
               continue;
            end
            th_a = asin(s_a/rho);
            z_a     = cos(th_a);
            z_min   = cos(pi);
            z_ar    = linspace(z_min, z_a, nz);

            rholam2= (rho/lam)^2    ;    

            pprr=-Iext*ri/(2*pi*d);

            z_ar    = linspace(z_min, z_a, nz);
            sol = pdepe(m,@(z,t,u,dudz)pdefun(z,t,u,dudz, tau,rholam2),@icfun,@(zL,uL,zR,uR,t)bcfun(zL,uL,zR,uR,t,pprr),z_ar,t_ar);
            va_vs_t = sol(:,nz); % v(t,th=theta_a)
            % for debug, only works with regular for, not parfor
            %figure(3); plot(t_ar,va_vs_t); hold on
            %figure(4); surf(sol)
            va_vs_t = va_vs_t /max(va_vs_t);
            it=find(va_vs_t>0.5*max(va_vs_t),1);
            if (it<10)
                printf("reaching time res, increase nt\n")
            end        
            t50(j,i) = t_ar(it); % 50% time
            vss(j,i) = max(sol(:,nz)); % input resistance
        end
    end
    

    figure(1)
        set(gcf,'position',[0,0,400,300])
        t50_exp= t50exp * ones(length(rho_ar),1);
        t50_erf= t50erf * ones(length(rho_ar),1);
        plot(rho_ar,t50','-','LineWidth',1)
        hold on
        plot(rho_ar,t50_erf','LineWidth',1);
        plot(rho_ar,t50_exp','LineWidth',1);

        legend(string(s_a_ar))
        xlabel('\rho')
        ylabel('Half maximum time (\tau)')
        ylim([0 0.75])
        annotation('textbox', [0.45, 0.4, 0.18, 0.05], 'String', 'Infinite Cable','EdgeColor','white')
        annotation('textbox', [0.45, 0.8, 0.4, 0.05], 'String', 'Single Compartment','EdgeColor','white')
        savefig('sphere_t50_rho_lam.fig')
        saveas(gcf,'sphere_t50_rho_lam.png')

    figure(2)    
        plot(rho_ar,vss','-','LineWidth',1)
        legend(string(s_a_ar))
         xlabel('\rho')
        ylabel('vss')

    figure(3) 
        vss_temp=vss.*(4*pi*rho_ar.^2)';
        plot(rho_ar,vss_temp','-','LineWidth',1)
        legend(string(s_a_ar))
        xlabel('\rho')
        ylabel('vss X surfaceArea')

        
    t50=t50';
    t50p=[s_a_ar' t50 t50_exp t50_erf];

    csvwrite('t50_sphere_rho_lam_case2.csv',t50p)

    csvwrite('vss_sphere_rho_lam_case2.csv',vss)

case(3) % Single run, square pulse,

    rho_ar = [lam/55];
    tha_ar = 0.025;

    nrho=length(rho_ar);
    nth=length(tha_ar);

    t50=zeros(nrho,nth);
    vss=zeros(nrho,nth);
    for j=1:nrho
        rho = rho_ar(j);
        rholam2= (rho/lam)^2        
        
        pprr=-Iext*ri/(2*pi*d);
        ton=0;
        toff=10*tau;
        
        for i=1:nth
            th_a = tha_ar(i);
            z_a     = cos(th_a);
            z_min   = cos(pi);
            z_ar    = linspace(z_min, z_a, nz);
            sol = pdepe(m,@(z,t,u,dudz)pdefun(z,t,u,dudz, tau,rholam2),@icfun,@(zL,uL,zR,uR,t)bcfun_pulse(zL,uL,zR,uR,t,ton,toff,pprr),z_ar,t_ar);
            va_vs_t = sol(:,nz); % v(t,th=theta_a)
            % for debug, only works with regular for, not parfor
            %figure(3); plot(t_ar,va_vs_t); hold on
            vss=sol(ceil(toff/dt-1),:);
            
            sol_plot=downsample(downsample(sol,ceil(nt/100))',ceil(nz/100))';
            figure(1)
                surf(sol_plot)
            
            figure(2)
                plot(t_ar,va_vs_t)
                ylabel('V at pipette')
                xlabel('time')
                hold on
                plot(t_ar_plot,exp_model*max(va_vs_t))
                
            figure(3)
                plot(vss);
                xlabel('position  as cos(theta)')
                ylabel('voltage at peak')
                
          %  va_vs_t = va_vs_t /max(va_vs_t);
            it=find(va_vs_t>0.5*max(va_vs_t),1);
            if (it<10)
                printf("reaching time res, increase nt\n")
            end        
            t50(j,i) = t_ar(it) % 50% time
            t50rel =  t_ar(it)/t50exp
            vss(j,i) = max(sol(:,nz)); % input resistance
        end
    end
    
    
end % switch
