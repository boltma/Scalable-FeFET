function [delP] = Switched_single(a, b, p, q, Pr, tauo, alpha, bet, amp, pw)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    tau = @(Ea) tauo * exp((Ea/amp).^alpha); %time constant
    sp = @(Ea) 1 - exp(-(pw./tau(Ea)).^bet); %switching probability
    f_Ea = @(Ea) (a/b) * (Ea/b).^(a*p-1) / beta(p,q) ./ (1 + (Ea/b).^a).^(p+q);
    
    sportion = @(Ea) f_Ea(Ea) .* sp(Ea);
    delP = -Pr + 2 * Pr * integral(sportion, 0, Inf);
%     delP = 2 * Pr * integral(sportion, 0, Inf);
end

