function [eigA, eigAn, eigV, eigVn] = func_classify_singular(G, H, Hn, w, mu_step, weight)
    
    if nargin <= 5
        weight = 0.3;
    end

    G = evalfr(G, 1i*w);

    n = size(H, 1)/2;

    Sn = inv(eye(2*n) + G * Hn);
    eigAn = eig(Sn(1:n, 1:n));
    eigVn = eig(Sn(n+1:end, n+1:end));

    mu_arr = linspace(0, 1, 1/mu_step);

    eig_arr = zeros(length(mu_arr), 2*n);

    for rr = 1:length(mu_arr)
        H_mu = (1 - mu_arr(rr)) * Hn + mu_arr(rr) * H;
        R_mu = eye(2*n) + G * H_mu;

        eig_arr(rr,:) = 1./eig(R_mu);
    end

    eig_arr_re = func_track_eigenvalues(eig_arr, weight);

    tol = 1e-5;
    eig_idx_A = [];
    eig_idx_V = [];
    for ii = 1:size(eig_arr_re, 2)
        err_A = abs(eigAn - eig_arr_re(1,ii))./abs(eig_arr_re(1,ii));
        err_V = abs(eigVn - eig_arr_re(1,ii))./abs(eig_arr_re(1,ii));

        if min(err_A) < min(err_V) && length(eig_idx_A)<n % length(eig_idx_A)<n to ensure length(eig_idx_A)==length(eig_idx_V)
            eig_idx_A = [eig_idx_A, ii];
        else
            eig_idx_V = [eig_idx_V, ii];
        end

        % index_A = find(err_A < tol);
        % index_V = find(err_V < tol);
        % if ~isempty(index_A)
        %     eig_idx_A = [eig_idx_A, ii];
        % elseif ~isempty(index_V)
        %     eig_idx_V = [eig_idx_V, ii];
        % end
    end

    eigAn = eig_arr_re(1, eig_idx_A);
    eigA  = eig_arr_re(end, eig_idx_A);

    eigVn = eig_arr_re(1, eig_idx_V);
    eigV  = eig_arr_re(end, eig_idx_V);

end