function [MW, Vth1, Vth2] = get_MW_simul(amp, pw, St_init, Weight, r_Ea, r_voff, Pr, tauo, alpha, bet, epife, Ndom)
% function [vfev,Pcurr] = get_MW(amp, pw, Pr, tauo, alpha, bet, epife)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

    %% first define the waveform
    [time,volt,rstend,stend] = wfdef_full(amp,pw,20e-9,1e-5,pw/300);
%     [time,volt,rstend,stend] = wfdef_full(amp,pw,max(pw/1000,10e-9),1e-5,pw/300);

    
%     %% finish the simulation
%     Ndom = 100; %number of domains
    
    % parameter
    T = 300; %temperature, K
    Na = 3e17; %substrate doping
    til = 1e-7; %interlayer thickness
    tfe = 0.8e-6; %ferroelectric thickness
    epiv = 8.85e-14;
    W = 1;
    L = 1;
    VFB = 0;
    tstep = 0;
    
    
    % start the simulation
    h = zeros(Ndom, 1); %intiialize the h
    hu = zeros(Ndom, 1); %intiialize the h
    tau = zeros(Ndom,1);
    vfe = 0;
    vfe_pre = 0;
    vfev = [];
    vgv = [];
    timv = [];
    Pcurr = [];
    Qcurr = [];
    h1test = [];
    tau1test = [];
    Stsum = [];
    
    % treat the voltage 0
    TIMELIMIT = 1e9;
    vswitchlimit = r_Ea/((log(TIMELIMIT/tauo))^(1/alpha));
    
    St = St_init;
    rng(0);
    for indii=1:length(time) %if the current time does not exceed the total time, keeps the simulation
        vcurr = volt(indii);
        vfe_pre = vfe;
        
        %% time step determination
        if indii == 1
            tsteptmp = time(2) - time(1);
        else
            tsteptmp = time(indii) - time(indii-1); 
        end
        
        tstep = tsteptmp;
        srand = rand(Ndom,1);
        
        F = @(x) Pstate(St,Weight, tstep, h, x, vfe_pre, r_voff, r_Ea, vswitchlimit, tauo, alpha, bet, srand) / sum(Weight) * Pr + 1e6 * epiv * epife / tfe * x - ...
            1e6 * MOSFET_Qmos(til, Na, T, W, L, vcurr-x, VFB, 0, 0); 
        vfe = fzero(F,vfe_pre);
        [St, hu] = Pstate_ret(St, tstep, h, vfe, vfe_pre, r_voff, r_Ea, vswitchlimit, tauo, alpha, bet, srand);

        vgv = [vgv volt(indii)];
        vfev = [vfev vfe];
        timv = [timv time(indii)];
        
        Stsum = [Stsum sum(St.*Weight)/sum(Weight)*Pr];
        
        h = hu;
        
    end
    MW = vfev(rstend) - vfev(stend);
    Vth1 = vfev(rstend);
    Vth2 = vfev(stend);
    

end

