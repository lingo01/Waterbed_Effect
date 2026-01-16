function eigenvalues_matrix = func_track_eigenvalues(eig_vals, weight)
    % 输入:
    % eig_vals: 维数为(num_params, num_eigenvalues) 的特征值矩阵

    % num_eigenvalues: 矩阵的阶数
    % num_params: 参数扫描的次数

    if nargin == 1
        weight = 0;
    end
    

    eig_vals = eig_vals';

    num_eigenvalues = size(eig_vals, 1);
    num_params = size(eig_vals, 2);
    
    eigenvalues_matrix = zeros(num_eigenvalues, num_params);

    % 初始匹配
    initial_eigenvalues = eig_vals(:, [1,2]);
    [~, idx] = sort(abs(initial_eigenvalues(:,1)));
    eigenvalues_matrix(:, 1) = initial_eigenvalues(idx,1);
    [~, idx] = sort(abs(initial_eigenvalues(:,2)));
    eigenvalues_matrix(:, 2) = initial_eigenvalues(idx,2);

   %  eigenvalues_matrix(:, [1,2]) = eigenvalues_matrix(:, [1,2]);

   if num_params >= 3
    % 对于后续参数p进行计算和匹配
        for k = 3:num_params
            current_eigenvalues = eig_vals(:, k);
            last_eigenvalues = eigenvalues_matrix(:, k-1);
            last2_eigenvalues = eigenvalues_matrix(:, k-2);
    
            % 计算距离矩阵
            distances_1 = abs(current_eigenvalues.' - last_eigenvalues);
            distances_2 = abs(current_eigenvalues.' - 2*last_eigenvalues + last2_eigenvalues);
    
            distances = weight * distances_1 + (1-weight) * distances_2;
            
            % 贪婪法匹配特征值：每次选择距离最小的匹配
            for i = 1:num_eigenvalues
                [~, idx] = min(distances(:));
                [row, col] = ind2sub(size(distances), idx);
                eigenvalues_matrix(row, k) = current_eigenvalues(col);
    
                % 将选中的行列设置为Inf，防止重复选择
                distances(row, :) = Inf;
                distances(:, col) = Inf;
            end
        end
   end

    eigenvalues_matrix = eigenvalues_matrix';
end