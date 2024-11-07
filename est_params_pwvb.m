function [lambda,t,prop_est,p,q] = est_params_pwvb(Azz, Azy, Ayy, Ayz, psi10, psi01, psi11)

    m=size(Azz,1);
    E=ones(m,m);
    e=ones(m,1);
    phi = psi10 + psi11;
    xi = psi01 + psi11;
    DM = E-eye(m);
    
    D_zy = diag(Azy);
    Azy = Azy - diag(diag(Azy));
    
    tp1 = (1-phi)'*Azz*(1-phi) + phi'*Azz*phi + (1-xi)'*Ayy*(1-xi) + xi'*Ayy*xi...
        +2*(1-phi)'*Azy*(1-xi) + 2*(phi)'*Azy*(xi)...  
        +2*(1-psi10-psi01)'*D_zy;
    tp2 = (1-phi)'*DM*(1-phi) + (phi)'*DM*(phi) + (1-xi)'*DM*(1-xi)...
        + (xi)'*DM*(xi) + 2*(1-phi)'*DM*(1-xi)+ 2*(phi)'*DM*(xi)...
        + 2*(1-psi10-psi01)'*e;
    p = tp1/tp2;
    
    tq1 = (1-phi)'*Azz*(phi) + phi'*Azz*(1-phi) + (1-xi)'*Ayy*(xi) + (xi)'*Ayy*(1-xi)...
        +2*(1-phi)'*Azy*(xi) + 2*(phi)'*Azy*(1-xi)...
        +2*(psi10+psi01)'*D_zy;
    tq2 = (1-phi)'*DM*(phi) + (phi)'*DM*(1-phi)+ (1-xi)'*DM*(xi)...
        + (xi)'*DM*(1-xi) + 2*(1-phi)'*DM*(xi) + 2*(phi)'*DM*(1-xi)...
        + 2*(psi10+psi01)'*e;
    q = tq1/tq2;
    
    
    t = 0.5*log(p*(1-q)/(q*(1-p)));
    lambda = 1/(2*t+1e-6) * log((1-q)/(1-p));
    prop_est=1/2;
end