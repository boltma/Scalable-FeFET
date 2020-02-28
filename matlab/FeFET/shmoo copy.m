Ndev = 5;
figureN = 0;
amp = 2:0.4:3;
pw = 1e-6:4e-6:11e-6;
Ids_on =  zeros(figureN, Ndev);
Ids_off = zeros(figureN, Ndev);
Vg_read = 1.0;
for i=1:length(amp)
    for j=1:length(pw)
        figureN = figureN+1;
        disp('figure #,  amp(V),  pw(us)');
        display([figureN, amp(i), pw(j)*1e6]);
        [Ids_on(figureN, :), Ids_off(figureN,:)] = FeFET_write_Id_func(Ndev, amp(i), pw(j), figureN, Vg_read);
    end
end

% k = find(min(IDset,[],2) > 1e-6); 
% % to ensure that read out current in the on-state is larger than 1uA
% 
% Ids_on = IDset(k(1),:);
% Ids_off = IDreset(k(1),:);