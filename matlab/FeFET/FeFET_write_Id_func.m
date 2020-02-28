function [Ids_on, Ids_off, Vg_read] =FeFET_write_Id_func(Ndev, amp, pw, figureN)
%% v4
%% this script is used to calculate the Id-Vg curves after program and erase
%Ids_on =  zeros(Ndev);
%Ids_off = zeros(Ndev);
%%
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
voffset = 0;
coeff_init = [a b p q Pr tauo alpha bet epife voffset];

%% device parameters
kb = 8.6173303e-5; % boltzman constant, eV/K
T = 300; %temperature, K
Na = 3e17; %substrate doping
til = 1e-7; %interlayer thickness
tfe = 0.8e-6; %ferroelectric thickness
epiv = 8.85e-14;
W = 1;
L = 1;
VFB = 0;
miu = 50;

%% setup
rng(0);
Ndom = 2000;
%Ndev = 100; %number of devices
%Ndev = 100; %number of devices
Weight1 = normrnd(1,0.27,[Ndom,1]);
Weight = Weight1.^2;
St_init = (randi([0 1], Ndom, Ndev)-0.5)*2;  

r_Ea = normrnd(a,b,[Ndom,Ndev]);
r_voff = normrnd(0, voffset, [Ndom,Ndev]);

% define the amplitude and pulse width for write
%amp = 2.5;
%amp = 2;
%pw = 1e-6;
delay = 1e-5;
tstep = 2e-8;
transit = 1e-9;

% define the waveform
[time,volt,resetind,setind]=wfdef_full(amp, pw, tstep, delay, transit);

%% measurement parameters
vgstart = -0.5;
vgend = 1.7;
vgstep = 0.02;
VG = [vgstart:vgstep:vgend];
VD = 0.05;
VS = 0;
VFB = 0;

%% simulate the Polarization states for Ndev devices
parfor ii=1:Ndev
%    fprintf('processing device %d\n', ii);
    [vfev(:,ii), Stsum(:,ii), h1test(:,ii), tau1test(:,ii)] = get_FE_state(time, volt,St_init(:,ii), Weight, r_Ea(:,ii), r_voff(:,ii), Pr,tauo,alpha, bet,  epife, Ndom);
end

parfor ii=1:Ndev
%    fprintf('measuring device %d\n', ii);
    %% then get the IdVg curves with the polarization state
    IDreset(:,ii) = get_ID(Stsum(resetind,ii), epife, tfe, til, miu, Na, T, W, L, VG, -0.5, VD, VS);
    IDset(:,ii) = get_ID(Stsum(setind,ii), epife, tfe, til, miu, Na, T, W, L, VG, -0.5, VD, VS);
end

% k = find(abs(VG-Vg_read) < 0.005);
% Ids_on  = IDset(k,:);
% Ids_off = IDreset(k,:);

k = find(min(IDset,[],2) > 1e-6); 
%To ensure that VG index for read out current in the on-state > 1uA
Ids_on = IDset(k(1),:);    %get on-state current
Ids_off = IDreset(k(1),:); %get off-state current at the minimum Vg_read
Vg_read = VG(k(1));        %get the gate read voltage

figure(figureN);
%semilogy(VG,IDreset,'r')
plot(VG,IDreset,'r');
hold on;
%semilogy(VG,IDset,'b')
plot(VG,IDset,'b');
title(strcat('amplitude = ',num2str(amp),'V, width = ', num2str(pw,'%e'), 'V, figure #: ', num2str(figureN)));
end