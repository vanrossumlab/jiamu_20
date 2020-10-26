function [c,f,s] = pdefun(z,t,u,dudz,tau,rholam2)    
    c = tau*rholam2;
    f = (1-z^2)*dudz;
    s = -rholam2*u;
end
