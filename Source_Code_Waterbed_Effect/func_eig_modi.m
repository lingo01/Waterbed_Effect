function eig_out = func_eig_modi(mode, mat)
    n = size(mat, 1);
    if strcmp(mode, 'matS')
        matL = inv(mat) - eye(n);
    elseif strcmp(mode, 'matR')
        matL = mat - eye(n);
    end

    eigR = eig(matL) + ones(n, 1);

    if strcmp(mode, 'matS')
        eig_out = 1./eigR;
    elseif strcmp(mode, 'matR')
        eig_out = eigR;
    end
end