
function [lambda,t,prop_est,p,q]=est_params(P,pi)
n=size(P,1);
E=ones(n,n);
e=ones(n,1);
p=(pi'*P*pi+(1-pi)'*P*(1-pi))/(pi'*(E-eye(n))*pi+(1-pi)'*(E-eye(n))*(1-pi));
q=((1-pi)'*P*pi)/((1-pi)'*(E-eye(n))*pi);


t = 0.5 * log(p*(1 - q)/(q*(1 - p))); lambda = log((1 - q)/(1 - p))/(2*t+1e-6);
prop_est=1/2;

end
