function [r_Ea] = Ea_dist(a, b, p, q, Ndom)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    Ea = linspace(0,8,Ndom);
    f_Ea = (a/b) * (Ea/b).^(a*p-1) / beta(p,q) ./ (1 + (Ea/b).^a).^(p+q);
    r_Ea = randpdf(f_Ea,Ea,[Ndom,1]);
    figure(1); hold on;
    histogram(r_Ea)
end

