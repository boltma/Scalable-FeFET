function [delP] = Switched_full(coeff,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    a = coeff(1);
    b = coeff(2);
    p = coeff(3);
    q = coeff(4);
    Pr = coeff(5);
    tauo = coeff(6);
    alpha = coeff(7);
    bet = coeff(8);

    %input amp and pw
    amp = x(:,2);
    pw = x(:,1);
    delP = zeros(length(pw),1);
    
    for ii=1:length(pw) %get the switched polarization for each pulse parameter
        delP(ii) = Switched_single(a, b, p, q, Pr, tauo, alpha, bet, amp(ii), pw(ii));
    end
end

