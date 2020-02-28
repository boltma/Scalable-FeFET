function [ Qmos ] = MOSFET_Qmos( tox, Na, temp, W, L, VGB, VFB, VDB, VSB )
%this function calculates the drain current of MOSFET
%this expression contains drift and diffusion current
%the reference is pp. 162, "Operation and Modeling of the MOS transistor"
    %% define some constants
    n_i = 5.29 * 1e19 * (temp / 300)^2.54 * exp(-6726 / temp); %unit cm-3
    k = 1.380648e-23;             % Boltzmann constant (J/K)
    qe = 1.6e-19;
    Vt = k * temp / qe;   %thermal voltage (V)
    
    epi_v = 8.85e-14; %F/cm
    epi_s = 11.8;
    epi_o = 3.9;
    
    LD = sqrt(Vt * epi_s * epi_v / Na / qe ); %cm
    Cox = epi_o * epi_v / tox; %F/cm2
    phib = Vt * log( Na /n_i); %V
    gamma = sqrt(2 * qe * Na * epi_s * epi_v) / Cox;  %V^1/2
    
    %% first solve for phis and phid, given the boundary conditions
    F_S = @(phi) sqrt( (exp(-1 * phi / Vt) + phi / Vt - 1 ) + (n_i/Na)^2 * exp(- VSB / Vt) * ( exp( phi / Vt) - phi / Vt - 1 - ((phi / Vt).^2 ./ ((phi / Vt).^2 + 2)) ) );
    F_D = @(phi) sqrt( (exp(-1 * phi / Vt) + phi / Vt - 1 ) + (n_i/Na)^2 * exp(- VDB / Vt) * ( exp( phi / Vt) - phi / Vt - 1 - ((phi / Vt).^2 ./ ((phi / Vt).^2 + 2)) ) );
    Eqn_S = @(phi) sqrt(2) * sign(phi) .* epi_s .* epi_v .* Vt ./ LD .* F_S(phi) / Cox + phi - (VGB - VFB);
    Eqn_D = @(phi) sqrt(2) * sign(phi) .* epi_s .* epi_v .* Vt ./ LD .* F_D(phi) / Cox + phi - (VGB - VFB);
    phis = real(fzero(Eqn_S, phib));
    phid = real(fzero(Eqn_D, phib));
    
%     %second solve the drain current
%     IDdrift = W / L * miu * Cox * ( (VGB - VFB) * (phid - phis) - 1/2 * (phid^2 - phis^2) - 2/3 * gamma * ( phid^(3/2) - phis^(3/2) ) );
%     IDdiff  = W / L * miu * Cox * ( Vt * (phid - phis) + Vt * gamma * ( phid^(1/2) - phis^(1/2) ));
%     ID = IDdrift + IDdiff;0.1
    
    %third solve the total gate charge
    %the reference is the paper "SP:an advanced surface potential based compact MOSFET model"
    phim = (phis + phid) / 2.0;
    dphi = phid - phis;
    if VGB < VFB || abs(VGB - VFB) <= 1e-9
        Qmos = Cox * W * L * (VGB - VFB - phim);
        qim = 0;
        alpha = 0;
    else
        alpha = 1 + gamma / 2 / sqrt(phim);
        qim = (VGB - VFB - phim - gamma * sqrt(phim));
        H = qim / alpha + Vt;
        Qmos = Cox * W * L * (VGB - VFB - phim + dphi^2/12/H);
    end
    Qmos = Qmos / W / L; %per unit area
end

