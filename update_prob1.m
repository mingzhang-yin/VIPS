function [psi10, psi01, psi11, psi00] = update_prob1(theta10,theta01,theta11)
    psi10 = 1./(1+exp(-theta10)+exp(theta01-theta10)+exp(theta11-theta10));
    psi01 = 1./(1+exp(-theta01)+exp(theta10-theta01)+exp(theta11-theta01));
    psi11 = 1./(1+exp(-theta11)+exp(theta10-theta11)+exp(theta01-theta11));
    psi00 = exp(-theta10)./(1+exp(-theta10)+exp(theta01-theta10)+exp(theta11-theta10));
end