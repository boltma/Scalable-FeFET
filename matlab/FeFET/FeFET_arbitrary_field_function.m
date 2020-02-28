function [Pcurr] = FeFET_arbitrary_field_function(coeff, input)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%% this script is used to model the ferroelectric film with the time domain NLS model
% it is mainly based on the model shown in paper "Monte Carlo Simulation of Switching Dynamics in Polycrystaline Ferroelectric Capacitors"
% it is used to simulate arbitrary waveform

    %% Initialize the parameters
    Ndom = 2000; %number of domains

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
    tfe = 1e-6;
    epiv = 8.85e-14; 
    
    T = 300; %temperature, K
    Na = 3e17; %substrate doping
    til = 0.8e-7; %interlayer thickness
    W = 1;
    L = 1;
    VFB = 0;
    

    rng(0);
    %% activation field distribution
    Ea = linspace(0,8,Ndom);
    f_Ea = (a/b) * (Ea/b).^(a*p-1) / beta(p,q) ./ (1 + (Ea/b).^a).^(p+q);
    r_Ea = randpdf(f_Ea,Ea,[Ndom,1]);
    r_Ea(2) = 1.8;
    r_voff = normrnd(0, voffset, [Ndom,1]);
    
%     r_Ea = normrnd(a,b,[Ndom,1]);

    %% start the Monte Carlo simulation
    St_init = (randi([0 1], Ndom, 1)-0.5)*2; %initilize the polarization direction, which is random
%     St_init = -1*ones(Ndom,1);
    St = St_init;
    
    %% input parameter
    time = input(:,1);
    volt = input(:,2);

    Pcurr = zeros(length(time),1);
    vpre = volt(1);

    h = zeros(Ndom, 1); %intiialize the h
    hu = zeros(Ndom, 1); %intiialize the h
    tau = zeros(Ndom,1);
    tautest1 = [];
    tautest2 = [];
    tautest3 = [];
    tautest4 = [];
    tautest5 = [];
 
    htest1 = [];
    htest2 = [];
    htest3 = [];
    htest4 = [];
    htest5 = [];
    
    Stsum = [];
    Stest1 = [];
    Stest2 = [];
    Stest3 = [];
    Stest4 = [];
    Stest5 = [];
    vfe = 0;
    vfe_pre = 0;
    vfev = [];
    

    TIMELIMIT = 1e9;
    vswitchlimit = r_Ea/((log(TIMELIMIT/tauo))^(1/alpha));

    for indii=1:length(time) %if the current time does not exceed the total time, keeps the simulation
        vcurr = volt(indii);
        F = @(x) sum(St) / Ndom * Pr + 1e6 * epiv * epife / tfe * x - 1e6 * MOSFET_Qmos(til, Na, T, W, L, vcurr-x, VFB, 0, 0); 
        vfe_pre = vfe;
        vfe = fzero(F,vfe_pre);
        vfev = [vfev vfe];
        %% time step determination
        if indii == 1
            tstep = time(2) - time(1);
        else
            tstep = time(indii) - time(indii-1); 
        end


        %% switching probability section
        vswitch = (vfe + vfe_pre)/2.0 - r_voff;
%         slope = (vcurr - vpre ) / tstep;
        taus = tauo * exp((r_Ea./max(abs(vswitch), vswitchlimit)).^alpha); %time constant for each domain

        tau = taus;

        hu(vswitch .* St <= 0) = h(vswitch .* St <= 0) + tstep ./ tau(vswitch .* St <= 0);
        hu(vswitch .* St > 0) = h(vswitch .* St > 0) - tstep ./ tau(vswitch .* St > 0);
        
        Pswi = 1 - exp(h.^bet - hu.^bet);
        Pswi(h>hu) = -0.1;
        hu(hu<0) = 0; %reset the h for those domains not switched

%         tautest1 = [tautest1 tau(1)];
%         tautest2 = [tautest2 tau(2)];
%         tautest3 = [tautest3 tau(3)];
%         tautest4 = [tautest4 tau(4)];
%         tautest5 = [tautest5 tau(5)];
%         
%         htest1 = [htest1 hu(1)];
%         htest2 = [htest2 hu(2)];
%         htest3 = [htest3 hu(3)];
%         htest4 = [htest4 hu(4)];
%         htest5 = [htest5 hu(5)];
%         %generate a random number and determine if the domain will be
%         %switched or not
%         srand = rand(Ndom,1);
%         hu(Pswi > srand) = 0; %reset h for those domains already switched
        hu(Pswi > 0.5) = 0; %reset h for those domains already switched
%         St(Pswi > srand) = -1* St(Pswi > srand); %for thos switching probablity greater than rand, switch it
        St(Pswi > 0.5) = -1* St(Pswi > 0.5);
        
%         Stest1 = [Stest1 St(1)];
%         Stest2 = [Stest2 St(2)];
%         Stest3 = [Stest3 St(3)];
%         Stest4 = [Stest4 St(4)];
%         Stest5 = [Stest5 St(5)];
        Pcurr(indii) = Pr * sum(St)/Ndom + 1e6 * vfe * epife * epiv / tfe;
        h = hu; %update the h
        Stsum = [Stsum sum(St)/Ndom*Pr];
    end
%     figure(2); hold on;
%     plot(time,vfev);
end

