%% this script is used to model the ferroelectric film with the time domain NLS model
% it is mainly based on the model shown in paper "Monte Carlo Simulation of Switching Dynamics in Polycrystaline Ferroelectric Capacitors"
clc;
clear;
%% Initialize the parameters
Ndom = 1000; %number of domains

% parameter set for data in literature
% %% activation field distribution
% a = 6.7582;
% b = 1.6686; %MV/cm
% p = 2.0998;
% q = 1.3641;
% 
% %% polarization parameters
% Pr = 22.7186e-6; %C/cm2
% tauo = 446.33e-9; %switching time at inifitie time
% alpha = 4.2408; % voltage acceleration factor
% bet = 1.9134; %switching probablity factor


%parameter set for our own data
%% activation field distribution
a = 6.9307;
b = 1.6698; %MV/cm
p = 0.5724;
q = 0.7239;

%% polarization parameters
Pr = 24.3608; %C/cm2
tauo = 3.5725e-6; %switching time at inifitie time
alpha = 3.2975; % voltage acceleration factor
bet = 2; %switching probablity factor

%% activation field distribution
Ea = linspace(0,4,Ndom);
f_Ea = (a/b) * (Ea/b).^(a*p-1) / beta(p,q) ./ (1 + (Ea/b).^a).^(p+q); %generalized beta distribution
r_Ea = randpdf(f_Ea,Ea,[Ndom,1]); % generated activation field random numbers

%% start the Monte Carlo simulation
% St = randi([0 1], Ndom, 1); %initilize the polarization direction, which is random
St_init = -1*ones(Ndom,1); % all the polarization are pointing downward

% define the input waveform
%experimental data
data = csvread('P_pw_our_own.csv',1,0);
pw = data(:,1);
amp = data(:,2);
delP_exp = data(:,3);

time = max(pw); %total 1ms
timestmp = logspace(-8, log10(time), 1000); %pulse width vector
timestep = diff([ 0 timestmp]);

% Vapp = linspace(0.6, 1.8, 13); %applied voltage
% Eapp = Vapp / 0.83; %field, unit MV/cm

Vapp = unique(amp); %unit V
Eapp = Vapp; %unit MV/cm, this is for thickness of 10nm

Pcurr = zeros(length(timestmp), length(Vapp)); %polarization matrix for different voltage and pulse width. unit uC/cm2
for ii=1:length(Vapp)
    St = St_init; %for each field, initialize the polarization at the begging
    tau = tauo * exp((r_Ea/Eapp(ii)).^alpha); %time constant for each domain
    h = zeros(Ndom, 1); %intiialize the h, which  tracks the history of polarization
    for jj=1:length(timestmp) %do the monte carlo simulation for each time step
        if jj == 1 %first time point
            hu = timestep(jj)./tau;
            h = zeros(Ndom, 1);
            Pswi = 1 - exp(h.^bet - hu.^bet);
        else
            hu = h + timestep(jj)./tau;
            Pswi = 1 - exp(h.^bet - hu.^bet);
        end
        
        %generate a random number and determine if the domain will be
        %switched or not
        srand = rand(Ndom,1);
        St(Pswi > srand) = 1; %for thos switching probablity greater than rand, switch it
        Pcurr(jj,ii) = Pr * sum(St)/Ndom;
        h = hu; %update the h
    end
end

%% plot the data

%plot
color_plot = ['k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c']; 
figure
num_amp = length(unique(amp));
num_pw = length(amp)/num_amp;
for ii = 1:num_amp
    Plot_P_experiment = delP_exp((ii-1)*num_pw+1:ii*num_pw);
    semilogx(pw(1:num_pw),Plot_P_experiment,'marker','o','linestyle','none','markersize',12,'color',color_plot(ii)); hold on
    semilogx(timestmp,Pcurr(:,ii),'color',color_plot(ii)); hold on
end