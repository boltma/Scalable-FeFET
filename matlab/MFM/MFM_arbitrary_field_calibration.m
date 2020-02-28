%% this script is used to model the ferroelectric film with the time domain NLS model
% it is mainly based on the model shown in paper "Monte Carlo Simulation of Switching Dynamics in Polycrystaline Ferroelectric Capacitors"
% clc;
% clear;
% 
%% activation field distribution
a = 8;
b = 2; %MV/cm
p = 0.6775;
q = 0.8115;

%% polarization parameters
Pr = 18; %C/cm2
tauo = 1e-7; %s
alpha = 3.2975;
bet = 2;
epife = 50;
vfb = 0.5; %this the built in field distribution

%% input parameter
data = csvread('MFM_scheme3.csv',1,0);
time = data(:,1);
volt = data(:,2);
P_exp = 1e6*data(:,3);
input = [time volt];

coeff_init = [a b p q Pr tauo alpha bet epife vfb];

delP_sim = MFM_arbitrary_field_function(coeff_init, input);
figure
plot(volt,P_exp,'o');
hold on;
plot(volt,delP_sim,'-','LineWidth',2);

% %% try the identical pulses
% pnum = 10;
% step = 5e-9;
% pw = 4e-7;
% delay = 1e-5;
% transit = 1e-7;
% tcurr = 0;
% amp = -1.5;
% wf = [];
% time = [];
% for ii=1:pnum
%     %delay
%     time = [time tcurr + linspace(0,floor(delay/step)-1, floor(delay/step))];
%     wf = [wf zeros(1,floor(delay/step))];
%     
%     tcurr = tcurr + floor(delay/step);
%     %transit
%     time = [time tcurr + linspace(0,floor(transit/step)-1, floor(transit/step))];
%     wf = [wf linspace(0,amp, floor(transit/step))];
%     tcurr = tcurr + floor(transit/step);
%     
%     %pulse
%     time = [time tcurr + linspace(0,floor(pw/step)-1, floor(pw/step))];
%     wf = [wf linspace(amp,amp, floor(pw/step))];
%     tcurr = tcurr + floor(pw/step);
%     
%     %transit
%     time = [time tcurr + linspace(0,floor(transit/step)-1, floor(transit/step))];
%     wf = [wf linspace(amp,0, floor(transit/step))];
%     tcurr = tcurr + floor(transit/step);
%     
%     %delay
%     time = [time tcurr + linspace(0,floor(delay/step)-1, floor(delay/step))];
%     wf = [wf zeros(1,floor(delay/step))];
%     
%     tcurr = tcurr + floor(delay/step);
% end
% 
% input = [step * time' wf'];
% delP_sim = MFM_arbitrary_field_function(coeff_init, input);



