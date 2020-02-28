function [time,volt,rstend,stend] = wfdef_full(amp, pw, step, delay, transit)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

    amp_init = 4; %initialization pulse
    pw_init = 1e-6; %initialization pusle width
    
    volt = [];
    time = [];
    tcurr = 0;
    
    %% one initialization pulses
    [timetmp, volttmp] = wfdef_single(amp_init,pw_init,pw_init/30,delay,transit);
    time = [time tcurr+timetmp];
    volt = [volt volttmp];
    tcurr = time(end);
    
    %% one reset pulse
    [timetmp, volttmp] = wfdef_single(-1*amp,pw,step,delay,transit);
    time = [time tcurr+timetmp(2:end)];
    volt = [volt volttmp(2:end)];
    tcurr = time(end);
    rstend = length(time);
    
    %% one set pulse
    [timetmp, volttmp] = wfdef_single(amp,pw,step,delay,transit);    
    time = [time tcurr+timetmp(2:end)];
    volt = [volt volttmp(2:end)];
    tcurr = time(end);
    stend = length(time);
    
%     %% one reset pulse
%     [timetmp, volttmp] = wfdef_single(-1*amp,pw,step,delay,transit);
%     time = [time tcurr+timetmp(2:end)];
%     volt = [volt volttmp(2:end)];
%     tcurr = time(end);
%     
%     %% one set pulse
%     [timetmp, volttmp] = wfdef_single(amp,pw,step,delay,transit);    
%     time = [time tcurr+timetmp(2:end)];
%     volt = [volt volttmp(2:end)];
%     tcurr = time(end);
end

