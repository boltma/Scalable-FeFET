function [MW, Vth1, Vth2] = MW_variation(coeff,input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %% activation field distribution
    a = coeff(1);
    b = coeff(2); %MV/cm
    p = coeff(3);
    q = coeff(4);

    %% polarization parameters
    Pr = coeff(5); %C/cm2
    tauo = coeff(6); %s
    alpha = coeff(7);
    bet = coeff(8);
    
    epife = coeff(9);
    voffset = coeff(10);
    Ndom = 2000;
    % activation field distribution
%     Ea = linspace(0,8,Ndom);
%     f_Ea = (a/b) * (Ea/b).^(a*p-1) / beta(p,q) ./ (1 + (Ea/b).^a).^(p+q);
%     r_Ea = randpdf(f_Ea,Ea,[Ndom,1]);

    %% initialization
    rng(0);
    Ndev = 100; %number of devices
    Weight1 = normrnd(1,0.27,[Ndom,1]);
    Weight = Weight1.^2;
    St_init = (randi([0 1], Ndom, Ndev)-0.5)*2;  

    r_Ea = normrnd(a,b,[Ndom,Ndev]);
    r_voff = normrnd(0, voffset, [Ndom,Ndev]);
    
    amp = input(:,2);
    pw = input(:,1);
    
    tic
%     index = [[7:12] [25:30]];
%     MW = zeros(1,length(pw));
%     for jj=1:length(index)
%     parfor ii=1:length(pw)
    for jj=1:Ndev
        fprintf('calculating device %d\n',jj);
        for ii=1:36
%             ii=index(jj)
%             ii
            [MW(ii,jj),Vth1(ii,jj),Vth2(ii,jj)] = get_MW_simul(amp(ii),pw(ii),St_init(:,jj), Weight, r_Ea(:,jj), r_voff(:,jj), Pr,tauo,alpha, bet,  epife, Ndom);
%             [MW(ii),Vth1(ii),Vth2(ii)] = get_MW_simul(amp(ii),pw(ii),St_init(:,1), r_Ea(:,1), r_voff(:,1), Pr,tauo,alpha, bet,  epife, Ndom);
        end
    end
    toc

    
    
end

