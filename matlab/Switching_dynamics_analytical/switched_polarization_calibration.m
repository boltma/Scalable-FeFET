%% this script is used to model the ferroelectric film with the time domain NLS model
% it is mainly based on the model shown in paper "Monte Carlo Simulation of Switching Dynamics in Polycrystaline Ferroelectric Capacitors"
clc;
clear;

%% activation field distribution
a = 12.1;
b = 1.79; %MV/cm
p = 0.691;
q = 0.633;

%% polarization parameters
Pr = 22.9; %C/cm2
tauo = 387e-9; %s
alpha = 4.11;
bet = 2.07;

%% input parameter
data = csvread('P_pw_our_own.csv',1,0);
% data = csvread('P_pw_literature.csv',1,0);
amp = data(:,2);
pw = data(:,1);
delP_exp = data(:,3);
input = [pw amp];

coeff_init = [a b p q Pr tauo alpha bet]; %initial guess of the parameters

lb=[1, 0.1, 0.1, 0.1, 10, 1e-9, 1, 2]; % lower bound constraints
ub=[100, 10  10, 10, 30, 1e-5, 10, 2]; % upper bound constraints
options = optimset('Display','iter','TolFun',1e-15, 'TolX', 1e-8);
for ii = 1:1
    c_pot = lsqcurvefit(@Switched_full, coeff_init, input, delP_exp,lb,ub,options) ;
    coeff_init = c_pot;
end %after this, c_pot gives the optimized parameter value

% plot the data
color_plot = ['k','r','b','m','g','c','k','r','b','m','g','c','k','r','b','m','g','c']; 
delP_sim = Switched_full(c_pot, input);
figure
num_amp = length(unique(amp));
num_pw = length(amp)/num_amp;
for ii = 1:num_amp
    Plot_P_experiment = delP_exp((ii-1)*num_pw+1:ii*num_pw);
    semilogx(pw(1:num_pw),Plot_P_experiment,'marker','o','linestyle','none','markersize',12,'color',color_plot(ii)); hold on
    Plot_P_sim = delP_sim((ii-1)*num_pw+1:ii*num_pw);
    semilogx(pw(1:num_pw),Plot_P_sim,'color',color_plot(ii)); hold on
end


