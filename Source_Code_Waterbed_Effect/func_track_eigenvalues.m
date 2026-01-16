function eigenvalues_matrix = func_track_eigenvalues(eig_vals, weight)
    % Inputs:
    % eig_vals: eigenvalue matrix of size (num_params, num_eigenvalues)
    %
    % num_eigenvalues: matrix order (number of eigenvalues per parameter point)
    % num_params: number of parameter samples

    if nargin == 1
        weight = 0;
    end
    

    eig_vals = eig_vals';

    num_eigenvalues = size(eig_vals, 1);
    num_params = size(eig_vals, 2);
    
    eigenvalues_matrix = zeros(num_eigenvalues, num_params);

    % Initial matching
    initial_eigenvalues = eig_vals(:, [1,2]);
    [~, idx] = sort(abs(initial_eigenvalues(:,1)));
    eigenvalues_matrix(:, 1) = initial_eigenvalues(idx,1);
    [~, idx] = sort(abs(initial_eigenvalues(:,2)));
    eigenvalues_matrix(:, 2) = initial_eigenvalues(idx,2);

   %  eigenvalues_matrix(:, [1,2]) = eigenvalues_matrix(:, [1,2]);

   if num_params >= 3
    % Match eigenvalues for subsequent parameter samples
        for k = 3:num_params
            current_eigenvalues = eig_vals(:, k);
            last_eigenvalues = eigenvalues_matrix(:, k-1);
            last2_eigenvalues = eigenvalues_matrix(:, k-2);
    
            % Compute distance matrices
            distances_1 = abs(current_eigenvalues.' - last_eigenvalues);
            distances_2 = abs(current_eigenvalues.' - 2*last_eigenvalues + last2_eigenvalues);
    
            distances = weight * distances_1 + (1-weight) * distances_2;
            
            % Greedy matching: repeatedly pick the smallest-distance pair
            for i = 1:num_eigenvalues
                [~, idx] = min(distances(:));
                [row, col] = ind2sub(size(distances), idx);
                eigenvalues_matrix(row, k) = current_eigenvalues(col);
    
                % Mask the selected row/column to avoid duplicate assignments
                distances(row, :) = Inf;
                distances(:, col) = Inf;
            end
        end
   end

    eigenvalues_matrix = eigenvalues_matrix';
end