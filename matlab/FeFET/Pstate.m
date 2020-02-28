function [sumSt] = Pstate(St_pre,Weight, tstep, h, vfe, vfe_pre, r_voff, r_Ea, vswitchlimit, tauo, alpha, bet, srand)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    St = St_pre;
    hu = h;
    vswitch = (vfe + vfe_pre)/2.0 - r_voff;
    tau = tauo * exp((r_Ea./max(abs(vswitch), vswitchlimit)).^alpha); %time constant for each domain
    hu(vswitch .* St <= 0) = h(vswitch .* St <= 0) + tstep ./ tau(vswitch .* St <= 0);
    hu(vswitch .* St > 0) = h(vswitch .* St > 0) - tstep ./ tau(vswitch .* St > 0);
    Pswi = 1 - exp(h.^bet - hu.^bet);
    Pswi(h>hu) = -0.1;
    hu(hu<0) = 0; %reset the h for those domains not switched
%     hu(Pswi > 0.5) = 0; %reset h for those domains already switched
%     St(Pswi > 0.5) = -1* St(Pswi > 0.5);
%     rng(0);
%     srand = rand(Ndom,1);
    hu(Pswi > srand) = 0;
    St(Pswi > srand) = -1* St(Pswi > srand);
    sumSt = sum(St.*Weight);
end

