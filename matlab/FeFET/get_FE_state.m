function [vfev, Stsum, h1test, tau1test] = get_FE_state(time, volt, St_init, Weight, r_Ea, r_voff, Pr, tauo, alpha, bet, epife, Ndom)
% function [vfev,Pcurr] = get_MW(amp, pw, Pr, tauo, alpha, bet, epife)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
    
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
%     rng(0);
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
        [St, hu, tau] = Pstate_ret(St,tstep, h, vfe, vfe_pre, r_voff, r_Ea, vswitchlimit, tauo, alpha, bet, srand);

        vgv = [vgv volt(indii)];
        vfev = [vfev vfe];
        timv = [timv time(indii)];
        
        Stsum = [Stsum sum(St.*Weight)/sum(Weight)*Pr];
        h1test = [h1test hu(14)];
        tau1test = [tau1test tau(14)];
        h = hu;
        
    end
    

end

