function [pL,qL,pR,qR] = bcfun(zL,uL,zR,uR,t,pprr)
    
    % BCs defined as : p+q*f(z,t,u,dudz)=0
    % here f=(1-z^2)*dudz= sin^2*dudz
    % V'(z_a) = -Iext *ri/(2*pi*d)/sin^2(za)

    pL = 0;
    
    qL = 1; % should be automatic
    
    pR = pprr; % -Iext*ri/(2*pi*d);
    
    qR = 1; 
end
