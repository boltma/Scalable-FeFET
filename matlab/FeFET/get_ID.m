function [ID] = get_ID(Stsum, epi_fe, tfe, til, miu, Na, T, W, L, VG, VFB, VD, VS)

    %get the drain current during the measurement
    Vfe = 0;
    Vfe_pre = 0;
    epiv = 8.85e-14;
    
    %get the length of VG
    for ii=1:length(VG)
        %solve the charge conservation equaiton
        F = @(x) Stsum/1e6 + epiv * epi_fe / tfe * x - MOSFET_Qmos(til, Na, T, W, L, VG(ii)-x, VFB, VD, VS); 
        Vfe_pre = Vfe;
        Vfe = fzero(F,Vfe_pre);
        ID(ii) = MOSFET_ID(til, miu, Na, T, W, L, VG(ii)-Vfe, VFB, VD, VS );
    end
        
end

