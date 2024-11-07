function [pi, P1, P2] = mf(P,pi_ini,S,lambda,t,prop,est_flag, v2)
n=size(P,1)/2;
half = 0.5*ones(2*n, 1);
N = size(P,1);

if est_flag==1
    [lambda,t,prop_est,p,q]=est_params(P,pi_old);
end

%initialize
pi_old = pi_ini;
prop_est = prop;
M = 4*t*(P - lambda*( ones(2*n) - eye(2*n) ) );
P1 = []; P2 = [];

for i = 1:S
    P1 = [P1, sum(pi_old)-N];
    P2 = [P2, pi_old' * v2];
    
    xi = M*(pi_old - half) + log(prop_est/(1-prop_est));      
    Pi(:, i) = 1./(exp(-xi) + 1);       
    pi_old = Pi(:, i);   

    
    if est_flag==1
        [lambda,t,prop_est,p,q]=est_params(P,pi_old);
        M = 4*t*(P - lambda*(ones(2*n) - eye(2*n)));
    end

end

pi=pi_old;
end
