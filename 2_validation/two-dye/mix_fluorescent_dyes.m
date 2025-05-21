function [mean_W1, mean_W2, std_W1, std_W2] = mix_fluorescent_dyes(mean_G_purelong, mean_S_purelong, ...
                                                                    mean_G_pureshort, mean_S_pureshort, ...
                                                                    G_values_filtered, S_values_filtered, numIndex)
    % Construct matrix A and vector B for solving weight equations
    A = [mean_G_purelong, mean_G_pureshort;
         mean_S_purelong, mean_S_pureshort;
         1, 1];
    
    % Get number of data points
    num_data_points = length(G_values_filtered);
    
    % Initialize matrix to store W1 and W2 results
    W_matrix = zeros(2, num_data_points);
    
    % Loop to solve for each weight vector w = [W1; W2]
    for i = 1:num_data_points
        B = [G_values_filtered(i); S_values_filtered(i); 1];
        W_matrix(:, i) = A \ B;
    end

    % Extract W1 and W2
    W1 = W_matrix(1, :);
    W2 = W_matrix(2, :);

    % Compute raw means for W1 and W2
    mean_W1 = mean(W1);
    mean_W2 = mean(W2);

    %% Thresholded normal distribution fitting for W1 and W2
    threshold = 0.4;  % 40% of peak height

    % --- W1 processing ---
    [W1_counts, W1_edges] = histcounts(W1, 30, 'Normalization', 'pdf');
    W1_centers = (W1_edges(1:end-1) + W1_edges(2:end)) / 2;

    [max_freq_W1, ~] = max(W1_counts);
    W1_threshold = max_freq_W1 * threshold;

    valid_indices_W1 = W1_counts >= W1_threshold;
    W1_filtered = W1_centers(valid_indices_W1);

    % Recompute mean and std from filtered values
    mean_W1 = mean(W1_filtered);
    std_W1 = std(W1_filtered);

    % Optional visualization (commented out)
    % histogram(W1, 15, 'Normalization', 'pdf'); hold on;
    % x_W1 = linspace(min(W1_filtered), max(W1_filtered), 100);
    % pdf_W1 = normpdf(x_W1, mean_W1, std_W1);
    % plot(x_W1, pdf_W1, 'r', 'LineWidth', 2); hold off;

    % --- W2 processing ---
    [W2_counts, W2_edges] = histcounts(W2, 30, 'Normalization', 'pdf');
    W2_centers = (W2_edges(1:end-1) + W2_edges(2:end)) / 2;

    [max_freq_W2, ~] = max(W2_counts);
    W2_threshold = max_freq_W2 * threshold;

    valid_indices_W2 = W2_counts >= W2_threshold;
    W2_filtered = W2_centers(valid_indices_W2);

    % Recompute mean and std from filtered values
    mean_W2 = mean(W2_filtered);
    std_W2 = std(W2_filtered);

    % Optional visualization (commented out)
    % histogram(W2, 15, 'Normalization', 'pdf'); hold on;
    % x_W2 = linspace(min(W2_filtered), max(W2_filtered), 100);
    % pdf_W2 = normpdf(x_W2, mean_W2, std_W2);
    % plot(x_W2, pdf_W2, 'r', 'LineWidth', 2); hold off;

end
