function [PF, D] = func_participant_factor(A)
    % PF(k,i)：第i个模态中第k个状态变量的参与程度
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