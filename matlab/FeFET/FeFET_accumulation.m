% this script is used for FeFET accumulation
clear;
clc;

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
Ndom = 20;

%% define the accumulation waveform
amp = 2.3;
pw = 1e-6;
[time,volt, index] = wfdef_acc(amp,pw,20e-9,1e-5,pw/300, 100);

rng(0);
Ndev = 1;
Ncycle = 10;

% Weight1 = normrnd(1,0.27,[Ndom,1]);
Weight1 = normrnd(1,0,[Ndom,Ndev]);
Weight = Weight1.^2;

St_init = (randi([0 1], Ndom, Ndev)-0.5)*2;  

r_Ea = normrnd(a,b,[Ndom,Ndev]);
r_voff = normrnd(0, vfb, [Ndom,Ndev]);

%% this is to generate device-to-device variation
% parfor ii=1:Ndev
%     fprintf('processing device %d\n', ii);
%     [vfev(:,ii), Stsum(:,ii), h1test(:,ii), tau1test(:,ii)] = get_FE_state(time, volt,St_init(:,ii), Weight, r_Ea(:,ii), r_voff(:,ii), Pr,tauo,alpha, bet,  epife, Ndom);
% end

%% for a single device, this is to generate cycle-to-cycle variation
jj = 1;
for ii=1:Ncycle
    tic
    fprintf('processing cycle %d\n', ii);
    [vfev(:,ii), Stsum(:,ii), h1test(:,ii), tau1test(:,ii)] = get_FE_state(time, volt,St_init(:,jj), Weight, r_Ea(:,jj), r_voff(:,jj), Pr,tauo,alpha, bet,  epife, Ndom);
    toc
end

%% find the flip of the polarization
for ii=1:Ncycle
    ind = find(Stsum(234:end, ii) > -1, 1);
    if isempty(ind)
        PN(ii) = nan;
    else
        PN(ii) = floor(ind/253);
    end
end

% %% calculate the IdVg curves
% % this is performed by getting the polarization state after each
% % accumulation pulse.
% 
% %once the polarization state is obtained, the Id-Vg is then simulated
% % device parameters
% kb = 8.6173303e-5; % boltzman constant, eV/K
% T = 300; %temperature, K
% Na = 3e17; %substrate doping
% til = 1e-7; %interlayer thickness
% tfe = 0.8e-6; %ferroelectric thickness
% epiv = 8.85e-14;
% W = 1;
% L = 1;
% VFB = 0;
% miu = 50;
% 
% % measurement parameters
% vgstart = -0.5;
% vgend = 1.7;
% vgstep = 0.02;
% VG = [vgstart:vgstep:vgend];
% VD = 0.05;
% VS = 0;
% VFB = 0;
% 
% for kk=1:Ncycle %50 IdVg curves
%     ID(kk,:) = get_ID(Stsum(index(kk),4), epife, tfe, til, miu, Na, T, W, L, VG, -0.5, VD, VS);
%     [uid, uidind, idind] = unique(ID(kk,:));
%     Vth(kk) = interp1(log10(uid+1e-20), VG(uidind), -7, 'linear', 'extrap');
% end

    