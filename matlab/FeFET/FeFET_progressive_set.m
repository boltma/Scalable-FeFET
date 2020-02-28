%FeFET_progressive_set
clear;
clc;
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
voffset = 0;

coeff_init = [a b p q Pr tauo alpha bet epife voffset];

Ndom = 20;
rng(0);
% activation field distribution
r_Ea = normrnd(a,b,[Ndom,1]);
r_voff = normrnd(0, voffset, [Ndom,1]);

%% initialization
% Weight1 = normrnd(1,0.27,[Ndom,1]);
Weight1 = normrnd(1,0,[Ndom,1]);
Weight = Weight1.^2;
St_init = (randi([0 1], Ndom, 1)-0.5)*2;

%% pulse parameters
ampstart = 2.85;
ampend = 2.85;
ampstep = 0.1;
pw = 1e-6;
delay = 1e-5;
tstep = 2e-8;
transit = 1e-9;

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

%% measurement parameters
vgstart = -0.5;
vgend = 1.7;
vgstep = 0.02;
VG = [vgstart:vgstep:vgend];
VD = 0.05;
VS = 0;
VFB = 0;
    
[time,volt, index] = wfdef_progressive_set(ampstart, ampend, ampstep, pw, tstep, delay, transit);
numcycle = 10;
tic
parfor ii=1:numcycle
    fprintf('processing cycle %d\n',ii);
    [vfev(:,ii),Stsum(:,ii)] = FeFET_simulation(time,volt, St_init, Weight, r_Ea, r_voff, Pr, tauo, alpha, bet, epife, Ndom); 
end
toc
index2 = index;
index2(length(index)+1) = length(time);
% use the polarization state after each pulse to get the ID-VG curves
for ii=1:numcycle
    for kk=1:length(index2)
        ID(kk,:,ii) = get_ID(Stsum(index2(kk),ii), epife, tfe, til, miu, Na, T, W, L, VG, -0.5, VD, VS);
        [uid, uidind, idind] = unique(ID(kk,:,ii));
        Vth(kk,ii) = interp1(log10(uid+1e-20), VG(uidind), -7, 'linear', 'extrap');
    end
end


    


