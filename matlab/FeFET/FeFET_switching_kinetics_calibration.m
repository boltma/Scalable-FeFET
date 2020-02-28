clc;
clear;

%% input data
data = csvread('MW_v2_clean_global_short.csv',1,0);
PW = data(:,1);    % Pulse width
Vapp = data(:,2);    % Pulse amplitude
MW_exp = data(:,3); % memory window

IN_X = [PW Vapp];

%% calibration parameters
% activation field distribution
% % a = 8;
% % b = 2.2; %MV/cm
% a = 2.2;
% b = 0.4; %MV/cm
% p = 0.6775;
% q = 0.8115;
% 
% 
% % polarization parameters
% Pr = 18; %C/cm2
% tauo = 4e-8; %s
% alpha = 2.5;
% bet = 4;
% epife = 22;
% vfb = 0;

% %% v2
% a = 2.2;
% b = 0.4; %MV/cm
% p = 0.6775;
% q = 0.8115;
% 
% 
% % polarization parameters
% Pr = 20; %C/cm2
% tauo = 2e-8; %s
% alpha = 2.5;
% bet = 3;
% epife = 25;
% vfb = 0;

%% v3
a = 2.3;
b = 0.4; %MV/cm
p = 0.6775;
q = 0.8115;


% polarization parameters
Pr = 25; %C/cm2
tauo = 1.9e-8; %s
alpha = 3.0;
bet = 2;
epife = 28;
vfb = 0;

coeff_init = [a b p q Pr tauo alpha bet epife vfb];

color_plot = ['k','r','b','m','g','c']; 
bias_num = 6;
PW_num = length(PW) / bias_num;
MW_fit = MW_variation(coeff_init,IN_X);
figure;
for ii = 1:bias_num
     pw_temp = PW( (ii-1)*PW_num + 1: (ii-1)*PW_num + PW_num );
     mw_temp = MW_fit( (ii-1)*PW_num + 1: (ii-1)*PW_num + PW_num,: );
     semilogx(pw_temp, mw_temp, 'color', color_plot(ii));
     hold on;
     mw_exp_temp = MW_exp( (ii-1)*PW_num + 1: (ii-1)*PW_num + PW_num );
     semilogx(pw_temp, mw_exp_temp, 'marker', 'o', 'linestyle','none','markersize',10, 'color', color_plot(ii));
     hold on;
end
xlim([1e-8, 1e-5]);

