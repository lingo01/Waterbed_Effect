function [idx_A, idx_V, idx_An, idx_Vn] = func_classify_track(eig_S, eig_Sn, eigA, eigAn, eigV, eigVn, w_row)
    eig_S  = eig_S(w_row, :);
    eig_Sn = eig_Sn(w_row, :);

    idx_A = zeros(size(eigA));
    idx_An = zeros(size(eigA));
    idx_V = zeros(size(eigA));
    idx_Vn = zeros(size(eigA));

    tol = 1e-3;

    for ii = 1:length(eigA)
        [~, idx_A(ii)] = min( abs(abs(eig_S)-abs(eigA(ii))) ./ abs(eigA(ii)) );
        [~, idx_An(ii)] = min( abs(abs(eig_Sn)-abs(eigAn(ii))) ./abs(eigAn(ii)) );
        [~, idx_V(ii)] = min( abs(abs(eig_S)-abs(eigV(ii))) ./abs(eigV(ii)) );
        [~, idx_Vn(ii)] = min( abs(abs(eig_Sn)-abs(eigVn(ii))) ./abs(eigVn(ii)) );
    end
end