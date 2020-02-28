function [time,volt, index] = wfdef_acc(amp, pw, step, delay, transit, cycle)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

    amp_init = -4; %initialization pulse
    pw_init = 1e-6; %initialization pusle width
    
    volt = [];
    time = [];
    tcurr = 0;
    
    %% one initialization pulses
    [timetmp, volttmp] = wfdef_single(amp_init,pw_init,pw_init/30,delay,transit);
    time = [time tcurr+timetmp];
    volt = [volt volttmp];
    tcurr = time(end);
    
    for ii=1:cycle  
        %% one reset pulse
        [timetmp, volttmp] = wfdef_single(amp+0*ii,pw,step,delay,transit);
        ind = find(volttmp ~= 0, 1) - 1;
        index(ii) = length(time) + ind -1;
        time = [time tcurr+timetmp(2:end)];
        volt = [volt volttmp(2:end)];
        tcurr = time(end);
    end

end

