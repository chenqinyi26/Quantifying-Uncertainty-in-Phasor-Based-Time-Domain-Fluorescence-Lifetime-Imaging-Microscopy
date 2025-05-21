function [G_values, G_values_filtered, S_values, S_values_filtered, ...
          Tau_values, Tau_values_filtered, ...
          total_photons_values, total_photons_values_filtered, ...
          histMatrix_filtered] = fluorescent_dyes(histAll)

    % Initialize new matrix
    global numSelectedBins StartBin
    [rows, cols, angles] = size(histAll);  % Get dimensions of input matrix
    
    numFiles = rows;
    numBins = cols;
    histMatrix_selected = zeros(numFiles, numSelectedBins);  % Preallocate new matrix
    
    % Iterate over each row and extract bins starting after StartBin
    for i = 1:numFiles
        startIndex = StartBin;  % Index after the intensity peak
        
        % Ensure we do not exceed matrix bounds
        endIndex = min(startIndex + numSelectedBins - 1, numBins);
        
        % Extract the desired segment into the selected matrix
        histMatrix_selected(i, :) = histAll(i, startIndex:endIndex);
    end
    
    % Load global timing variables
    global dt t_end f omega t
    
    %% Compute G, S, Tau for all pixels

    G_values = zeros(numFiles, 1);
    S_values = zeros(numFiles, 1);
    Tau_values = zeros(numFiles, 1);
    total_photons_values = zeros(numFiles, 1);
    
    for i = 1:numFiles
        % Compute total photon count
        total_photons = sum(histMatrix_selected(i,:));
        
        % Compute numerators for G and S
        G_numerator = sum(histMatrix_selected(i,:) .* cos(omega * t));
        S_numerator = sum(histMatrix_selected(i,:) .* sin(omega * t));
    
        if total_photons > 0
            G = G_numerator / total_photons;
            S = S_numerator / total_photons;
        else
            G = NaN;
            S = NaN;
            warning('Photon count is zero in file %s, G and S set to NaN', filename);
        end
        
        Tau = S / (omega * G);

        % Store results
        G_values(i) = G;
        S_values(i) = S;
        Tau_values(i) = Tau;
        total_photons_values(i) = total_photons;
    end

%% Outlier removal based on photon count histogram

threshold = 0.4; % Threshold for filtering based on histogram peak

% Compute histogram of total photon values
[total_photons_counts, total_photons_edges] = histcounts(total_photons_values, 'BinWidth', 10);
total_photons_centers = (total_photons_edges(1:end-1) + total_photons_edges(2:end)) / 2;

% Find the bin with highest frequency
[max_freq_total_photons, max_idx_total_photons] = max(total_photons_counts);
total_photons_peak = total_photons_centers(max_idx_total_photons);

% Determine threshold region
total_photons_threshold = max_freq_total_photons * threshold;
total_photons_valid_indices = find(total_photons_counts >= total_photons_threshold);
total_photons_min = total_photons_centers(min(total_photons_valid_indices));
total_photons_max = total_photons_centers(max(total_photons_valid_indices));

% Identify indices within valid range
filtered_indices = (total_photons_values >= total_photons_min & ...
                    total_photons_values <= total_photons_max);

% Extract filtered data
G_values_filtered = G_values(filtered_indices);
S_values_filtered = S_values(filtered_indices);
Tau_values_filtered = Tau_values(filtered_indices);
total_photons_values_filtered = total_photons_values(filtered_indices);

% Extract corresponding histograms
histMatrix_filtered = histMatrix_selected(filtered_indices, :);

% Display how many data points remain after filtering
disp(['Filtered data count: ', num2str(length(G_values_filtered))]);

end
