function [time,volt, index] = wfdef_progressive_set(ampstart, ampend, ampstep, pw, step, delay, transit)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

    amp_init = 6; %initialization pulse
    pw_init = 1e-6; %initialization pusle width
    
    volt = [];
    time = [];
    tcurr = 0;
    
    %%number of pulses
    Num = (ampend - ampstart)/ampstep + 1;
    
    [timetmp, volttmp] = wfdef_single(-1*amp_init,pw_init,pw_init/50,delay,transit);
    time = [time tcurr+timetmp];
    volt = [volt volttmp];
    tcurr = time(end);
    
    for ii=1:Num
        %% one initialization pulses
%         [timetmp, volttmp] = wfdef_single(-1*amp_init,pw_init,pw_init/50,delay,transit);
%         time = [time tcurr+timetmp];
%         volt = [volt volttmp];
%         tcurr = time(end);

        %% one reset pulse
        [timetmp, volttmp] = wfdef_single(ampstart+(ii-1)*ampstep,pw,step,delay,transit);
        ind = find(volttmp ~= 0, 1) - 1;
        index(ii) = length(time) + ind -1;
        time = [time tcurr+timetmp(2:end)];
        volt = [volt volttmp(2:end)];
        tcurr = time(end);
    end
    
end

