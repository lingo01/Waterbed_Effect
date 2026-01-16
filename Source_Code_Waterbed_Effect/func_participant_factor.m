function [PF, D] = func_participant_factor(A)
    % PF(k,i): participation factor of state k in mode i
    [V,D] = eig(A);
    D = diag(D);

    Phi = V; 
    Psi = inv(V);

    PF = zeros(size(V));
    for kk = 1:size(V, 1)
        for ii = 1:size(V, 2)
            PF(kk, ii) = Phi(kk, ii) * Psi(ii, kk);
        end
    end

end