function [time,volt] = wfdef_single(amp, pw, tstep, delay, transit)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

%% define the pulse
    volt = [];
    time = [];
    tcurr = 0;
    EPIS = 1e-15;
%     N = max(floor(1e-6/pw),1);
%     N = 33e-9/tstep;
    N = 100e-9/tstep;
    
    time = [time tcurr + N*tstep*(linspace(0,floor(delay/N/tstep), floor(delay/N/tstep)+1))];
    volt = [volt zeros(1,floor(delay/N/tstep)+1)];
    
    if abs(N*tstep*(floor(delay/N/tstep)) - delay) > EPIS
        time = [time delay];
        volt = [volt 0];
    end

    tcurr = tcurr + delay;
    %transit
    time = [time tcurr + tstep*(linspace(1,floor(transit/tstep), floor(transit/tstep)))];
    volt = [volt linspace(0,amp, floor(transit/tstep))];
    tcurr = tcurr + transit;
    
    if abs(tstep*(floor(transit/tstep)) - transit) > EPIS
        time = [time tcurr];
        volt = [volt amp];
    end

    
    %pulse
    time = [time tcurr + tstep * (linspace(1,floor(pw/tstep), floor(pw/tstep)))];
    volt = [volt linspace(amp,amp, floor(pw/tstep))];
%     tcurr = tcurr + floor(pw/tstep);
    tcurr = tcurr + pw;
    
    if abs(tstep*(floor(pw/tstep)) - pw) > EPIS
        time = [time tcurr];
        volt = [volt amp];
    end
    
    %transit
    time = [time tcurr + tstep*(linspace(1,floor(transit/tstep), floor(transit/tstep)))];
    volt = [volt linspace(amp,0, floor(transit/tstep))];
%     tcurr = tcurr + floor(transit/tstep);
    tcurr = tcurr + transit;
    
    if abs(tstep*(floor(transit/tstep)) - transit) > EPIS
        time = [time tcurr];
        volt = [volt 0];
    end
    
    %delay
    time = [time tcurr + N*tstep*(linspace(1,floor(delay/N/tstep), floor(delay/N/tstep)))];
    volt = [volt zeros(1,floor(delay/N/tstep))];
    tcurr = tcurr + delay;
    if abs(N*tstep*(floor(delay/N/tstep)) - delay) > EPIS
        time = [time tcurr];
        volt = [volt 0];
    end
    
end

