function [Pcurr] = MFM_arbitrary_field_function(coeff, input)
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

    rng(0);
    %% activation field distribution
    Ea = linspace(0,8,Ndom);
    f_Ea = (a/b) * (Ea/b).^(a*p-1) / beta(p,q) ./ (1 + (Ea/b).^a).^(p+q);
    r_Ea = randpdf(f_Ea,Ea,[Ndom,1]);
    r_voff = normrnd(0, voffset, [Ndom,1]);
    
%     r_Ea = normrnd(a,b,[Ndom,1]);

    %% start the Monte Carlo simulation
    % St = randi([0 1], Ndom, 1); %initilize the polarization direction, which is random
    St_init = ones(Ndom,1);
    St = St_init;
    
    %% input parameter
    time = input(:,1);
    volt = input(:,2);

    Pcurr = zeros(length(time),1);
    vpre = volt(1);

    h = zeros(Ndom, 1); %intiialize the h
    hu = zeros(Ndom, 1); %intiialize the h
    tau = zeros(Ndom,1);
    tautest = [];
    htest = [];
    Stsum = [];
    Stest1 = [];
    Stest2 = [];
    Stest3 = [];
    Stest4 = [];
    Stest5 = [];
    

    TIMELIMIT = 1e9;
    vswitchlimit = r_Ea/((log(TIMELIMIT/tauo))^(1/alpha));

    for indii=1:length(time) %if the current time does not exceed the total time, keeps the simulation
        vcurr = volt(indii);
        %% time step determination
        if indii == 1
            tstep = time(2) - time(1);
        else
            tstep = time(indii) - time(indii-1); 
        end


        %% switching probability section
        vswitch = (vcurr + vpre)/2.0 - r_voff;
%         slope = (vcurr - vpre ) / tstep;
        taus = tauo * exp((r_Ea./max(abs(vswitch), vswitchlimit)).^alpha); %time constant for each domain

        tau = taus;

        hu(vswitch .* St <= 0) = h(vswitch .* St <= 0) + tstep ./ tau(vswitch .* St <= 0);
        hu(vswitch .* St > 0) = h(vswitch .* St > 0) - tstep ./ tau(vswitch .* St > 0);
        
        Pswi = 1 - exp(h.^bet - hu.^bet);
        Pswi(h>hu) = -0.1;
        hu(hu<0) = 0; %reset the h for those domains not switched

        tautest = [tautest tau(123)];
        htest = [htest hu(123)];
        %generate a random number and determine if the domain will be
        %switched or not
        srand = rand(Ndom,1);
%         hu(Pswi > srand) = 0; %reset h for those domains already switched
        hu(Pswi > 0.5) = 0; %reset h for those domains already switched
        St(Pswi > srand) = -1* St(Pswi > srand); %for thos switching probablity greater than rand, switch it
        Stest1 = [Stest1 St(123)];
        Stest2 = [Stest2 St(475)];
        Stest3 = [Stest3 St(816)];
        Stest4 = [Stest4 St(881)];
        Stest5 = [Stest5 St(1)];
        Pcurr(indii) = Pr * sum(St)/Ndom + 1e6 * vcurr * epife * epiv / tfe;
        h = hu; %update the h
        vpre = vcurr;
        Stsum = [Stsum sum(St)/Ndom];
    end
end

