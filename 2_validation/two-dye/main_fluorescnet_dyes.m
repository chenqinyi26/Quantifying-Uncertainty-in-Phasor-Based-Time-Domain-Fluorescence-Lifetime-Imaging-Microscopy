clear all;
clc;

%% Load data
% Define the file path and filename template
filePath = 'histAll';
fileTemplate = 'histAll%d.mat';
filenum = 13;

% Define global variables
global  dt t_end f omega t numSelectedBins StartBin

dt = 9.77645305514160e-11;  % Time bin width
numSelectedBins = 40;       % Number of selected time bins
StartBin = 5;               % Start from a few bins after the max bin
t_end = dt * numSelectedBins;  % Total acquisition time
f = 1 / t_end;                 % Frequency set to 40 MHz
omega = 2 * pi * f;           % Angular frequency
t = dt/2:dt:t_end - (dt/2);

% Initialize cell arrays
G_values_all = cell(1, filenum);
S_values_all = cell(1, filenum);
Tau_values_all = cell(1, filenum);
G_values_filtered_all = cell(1, filenum);
S_values_filtered_all = cell(1, filenum);
Tau_values_filtered_all = cell(1, filenum);
total_photons_values_all = cell(1, filenum);
total_photons_values_filtered_all = cell(1, filenum);
histMatrix_filtered_all = cell(1, filenum);

% Loop through files 0 to filenum-1
for i = 0:filenum-1
    % Generate filename
    fileName = sprintf(fileTemplate, i);
    fullPath = fullfile(filePath, fileName);

    % Load the entire MAT file as a struct
    histData = load(fullPath);

    % Extract histAll_cleaned variable
    histAll = histData.histAll_cleaned; 

    % Call subfunction
    [G_on_purelong, G_values_filtered, ...
     S_on_purelong, S_values_filtered, ...
     Tau_values, Tau_values_filtered, ...
     total_photons_values, total_photons_values_filtered, ...
     histMatrix_filtered] = fluorescent_dyes_CNSI(histAll);

    % Store data in cell arrays
    G_values_all{i + 1} = G_on_purelong;
    S_values_all{i + 1} = S_on_purelong;
    Tau_values_all{i + 1} = Tau_values;
    total_photons_values_all{i + 1} = total_photons_values;
    G_values_filtered_all{i + 1} = G_values_filtered;
    S_values_filtered_all{i + 1} = S_values_filtered;
    Tau_values_filtered_all{i + 1} = Tau_values_filtered;
    total_photons_values_filtered_all{i + 1} = total_photons_values_filtered;
    histMatrix_filtered_all{i + 1} = histMatrix_filtered;
end

%% Pure components
for numIndex = [0, filenum - 1]  % Specify suffix indices for pure samples

    % Extract data from cell arrays
    G_values_filtered = G_values_filtered_all{numIndex + 1}; 
    S_values_filtered = S_values_filtered_all{numIndex + 1}; 
    Tau_values_filtered = Tau_values_filtered_all{numIndex + 1};
    histMatrix_filtered = histMatrix_filtered_all{numIndex + 1};

    % Generate candidate points on the universal semicircle
    theta = linspace(0, pi, 1000); % Angle range [0, pi]
    G_candidates = 0.5 + 0.5 * cos(theta);
    S_candidates = 0.5 * sin(theta);
    
    % Exclude (0,0) and (1,0)
    valid_idx = ~( (G_candidates == 0 & S_candidates == 0) | (G_candidates == 1 & S_candidates == 0) );
    G_candidates = G_candidates(valid_idx);
    S_candidates = S_candidates(valid_idx);
    
    % Compute total distances from candidates to all filtered G/S points
    total_distances = zeros(size(G_candidates));
    for i = 1:length(G_candidates)
        total_distances(i) = sum(sqrt((G_candidates(i) - G_values_filtered).^2 + (S_candidates(i) - S_values_filtered).^2));
    end
    
    % Find optimal point minimizing total distance
    [~, min_idx] = min(total_distances);
    G_on = G_candidates(min_idx);
    S_on = S_candidates(min_idx);
    
    fprintf('Optimal point: G = %.6f, S = %.6f\n', G_on, S_on);
    
    if numIndex == 0
        G_on_long = G_on;
        S_on_long = S_on;
    end
    if numIndex == filenum - 1
        G_on_short = G_on;
        S_on_short = S_on;
    end
end
Tau_on_long = S_on_long / (omega * G_on_long);
Tau_on_short = S_on_short / (omega * G_on_short);

fprintf('Tau_on_long = %.6f, G_on_long = %.6f, S_on_long = %.6f\n', Tau_on_long*1e9, G_on_long, S_on_long);
fprintf('Tau_on_short = %.6f, G_on_short = %.6f, S_on_short = %.6f\n', Tau_on_short*1e9, G_on_short, S_on_short);

%% Mixture weights
% Initialize result matrix
mix_W_matrix = zeros(filenum-2, 4); % Each row stores mean_W1, std_W1, mean_W2, std_W2

% Preallocate result storage
fitResults = cell(filenum-2, 2); % Store G and S fitting results

% Loop over dataset indices
for numIndex = 1:filenum-2
    % Extract data
    G_values_filtered = G_values_filtered_all{numIndex + 1}; 
    S_values_filtered = S_values_filtered_all{numIndex + 1}; 
    Tau_values_filtered = Tau_values_filtered_all{numIndex + 1};
    
    % Use phasor positions from the semicircle
    [mean_W1, mean_W2, std_W1, std_W2] = mix_fluorescent_dyes(G_on_long, S_on_long, ...
                                                              G_on_short, S_on_short, ...
                                                              G_values_filtered, S_values_filtered,...
                                                              numIndex);
    % Save results
    mix_W_matrix(numIndex, :) = [mean_W1, std_W1, mean_W2, std_W2];
end

%% Mixture uncertainty (from pure components)

% Initialize result arrays
Cov_uncertainty_matrix = cell(filenum-2, 1); 
G_uncertainty_array = zeros(filenum-2, 1); 
S_uncertainty_array = zeros(filenum-2, 1);

% Loop over numIndex = 1:9
for numIndex = 1:filenum-2
    % Compute mean photon count
    data_matrix = total_photons_values_filtered_all{numIndex + 1};
    total_photons = mean(data_matrix(:));

    % Extract weights
    P1 = mix_W_matrix(numIndex, 1);
    P2 = mix_W_matrix(numIndex, 3);

    % Use pure points from semicircle
    [Cov_matrix, G, S] = uncertainty_fluorescent_dyes(total_photons, P1, P2, ...
                                                      Tau_on_long, Tau_on_short);

    Cov_uncertainty_matrix{numIndex} = Cov_matrix;
    G_uncertainty_array(numIndex) = G;
    S_uncertainty_array(numIndex) = S;
end
%% Propagate G/S uncertainty to weight uncertainty

% Preallocate arrays to store weight distribution parameters
mu_W_all = zeros(filenum-2, 2);       % Mean weights for 9 datasets
Sigma_W_all = cell(filenum-2, 1);     % Covariance matrices for weights

% Pure component matrix (G, S) coordinates
pure_matrix = [G_on_long, G_on_short; S_on_long, S_on_short];

% Compute weight transformation matrix A
A = inv(pure_matrix); % Matrix for converting G-S coordinates to weights

% Iterate through 9 datasets
for i = 1:filenum-2
    % Extract current G-S mean
    mu_GS = [G_uncertainty_array(i); S_uncertainty_array(i)];

    % Extract current covariance matrix (2x2)
    Sigma_GS = Cov_uncertainty_matrix{i};

    % Compute weight mean and covariance
    mu_W = A * mu_GS;           % Weight mean
    Sigma_W = A * Sigma_GS * A'; % Weight covariance matrix

    % Store results
    mu_W_all(i, :) = mu_W';
    Sigma_W_all{i} = Sigma_W;
end

%% Merge data and plot 2D histogram

% Initialize merged arrays
G_combined = [];
S_combined = [];

% Indices to merge
indices_to_combine = [0,3,6,7,10,11,12];

% Merge selected datasets
for m = indices_to_combine
    G_combined = [G_combined; G_values_filtered_all{m + 1}];
    S_combined = [S_combined; S_values_filtered_all{m + 1}];
end

% Plot 2D histogram
numBins = 50;
[histCounts, edgesG, edgesS] = histcounts2(G_combined, S_combined, numBins, 'Normalization', 'count');

% Filter out low count bins
histCounts(histCounts < 1) = 0;

% Plot 2D histogram
figure;
h = histogram2('XBinEdges', edgesG, 'YBinEdges', edgesS, 'BinCounts', histCounts, ...
               'DisplayStyle', 'tile', 'Normalization', 'count');

xlabel('G', 'Color', 'k');
ylabel('S', 'Color', 'k');
title('Mixed Phasor Plot', 'Color', 'k');
colorbar;
axis equal;
xticks(0:0.1:1);
yticks(0:0.1:0.5);
caxis([1,20]);
set(gca, 'Box', 'on','XColor', 'k', 'YColor', 'k');
grid on;

% Plot semicircle
hold on;
theta = linspace(pi/2, pi, 100);
x = 0.5 + 0.5 * cos(theta);
y = 0.5 * sin(theta);
plot(x, y, 'black--', 'LineWidth', 1.5);
axis equal;

plot(G_on_long, S_on_long, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot(G_on_short, S_on_short, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
plot([G_on_long, G_on_short], [S_on_long, S_on_short], 'k--', 'LineWidth', 1.5);

% Draw uncertainty ellipses
ellipse_indices = [3,6,7,10,11];

hold on;
for i = ellipse_indices
    mean_GS = [G_uncertainty_array(i), S_uncertainty_array(i)];
    cov_GS = Cov_uncertainty_matrix{i};
    plot_gaussian_contour(mean_GS, cov_GS, 3, 'r-');
end
hold off;

%% G/S Uncertainty Distribution

figure('Position', [100, 100, 1200, 800]);
colors = lines(5);

% G distribution
subplot(2,1,1);
hold on;
G_all_values = [];
for idx = 1:length(ellipse_indices)
    numIndex = ellipse_indices(idx);
    G_all_values = [G_all_values; G_values_filtered_all{numIndex + 1}];
end
[G_counts, G_edges] = histcounts(G_all_values, 50);

exp_handles = gobjects(1, length(ellipse_indices));
uncertainty_handles = gobjects(1, length(ellipse_indices));

for idx = 1:length(ellipse_indices)
    numIndex = ellipse_indices(idx);
    G_values_filtered = G_values_filtered_all{numIndex + 1}; 
    G_uncertainty = G_uncertainty_array(numIndex);
    Cov_matrix = Cov_uncertainty_matrix{numIndex};
    sigma_G = sqrt(Cov_matrix(1,1));
    G_range = linspace(min(G_edges), max(G_edges), 500);
    G_pdf = normpdf(G_range, G_uncertainty, sigma_G);

    exp_handles(idx) = histogram(G_values_filtered, G_edges, 'Normalization', 'count', ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', colors(idx, :));

    uncertainty_handles(idx) = plot(G_range, G_pdf * max(G_counts) / max(G_pdf), '--', ...
        'Color', colors(idx, :) * 0.8, 'LineWidth', 2);
end

title('G Noisy Distribution with Uncertainty Model', 'Color', 'k');
xlabel('G', 'Color', 'k');
ylabel('Count', 'Color', 'k');
set(gca, 'Box', 'on','XColor', 'k', 'YColor', 'k');
grid on;
hold off;

% S distribution
subplot(2,1,2);
hold on;
S_all_values = [];
for idx = 1:length(ellipse_indices)
    numIndex = ellipse_indices(idx);
    S_all_values = [S_all_values; S_values_filtered_all{numIndex + 1}];
end
[S_counts, S_edges] = histcounts(S_all_values, 50);

for idx = 1:length(ellipse_indices)
    numIndex = ellipse_indices(idx);
    S_values_filtered = S_values_filtered_all{numIndex + 1}; 
    S_uncertainty = S_uncertainty_array(numIndex);
    Cov_matrix = Cov_uncertainty_matrix{numIndex};
    sigma_S = sqrt(Cov_matrix(2,2));
    S_range = linspace(min(S_edges), max(S_edges), 500);
    S_pdf = normpdf(S_range, S_uncertainty, sigma_S);

    histogram(S_values_filtered, S_edges, 'Normalization', 'count', ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', colors(idx, :));

    plot(S_range, S_pdf * max(S_counts) / max(S_pdf), '--', ...
        'Color', colors(idx, :) * 0.8, 'LineWidth', 2);
end

title('S Noisy Distribution with Uncertainty Model', 'Color', 'k');
xlabel('S', 'Color', 'k');
ylabel('Count', 'Color', 'k');
set(gca, 'Box', 'on','XColor', 'k', 'YColor', 'k');
grid on;
hold off;
%% Percentage of data points within 1σ, 2σ, 3σ ellipses

% Preallocate matrix: rows = number of selected ellipses, columns = σ ranges
percentage_matrix = zeros(length(ellipse_indices), 3);

for idx = 1:length(ellipse_indices)
    numIndex = ellipse_indices(idx);
    
    % Get data
    G_values_filtered = G_values_filtered_all{numIndex + 1};
    S_values_filtered = S_values_filtered_all{numIndex + 1};
    G_uncertainty = G_uncertainty_array(numIndex);
    S_uncertainty = S_uncertainty_array(numIndex);
    Cov_matrix = Cov_uncertainty_matrix{numIndex};
    
    % Calculate Mahalanobis distance
    data_points = [G_values_filtered - G_uncertainty, S_values_filtered - S_uncertainty];
    Cov_inv = inv(Cov_matrix);
    
    % Compute Mahalanobis squared distance d²
    mahalanobis_distances = sum((data_points * Cov_inv) .* data_points, 2);
    
    % Compute percentage of points within 1σ, 2σ, 3σ bounds
    percentage_matrix(idx, 1) = sum(mahalanobis_distances <= 1) / length(mahalanobis_distances) * 100; % 1σ
    percentage_matrix(idx, 2) = sum(mahalanobis_distances <= 4) / length(mahalanobis_distances) * 100; % 2σ
    percentage_matrix(idx, 3) = sum(mahalanobis_distances <= 9) / length(mahalanobis_distances) * 100; % 3σ
end

% Display results
disp('Percentage of data points within 1σ, 2σ, 3σ ellipses (unit: %):');
disp(percentage_matrix);

%% Compute weights for each data point and plot histograms

figure('Position', [100, 100, 2000, 250]); % Canvas size
tiledlayout(1,2); % 1 row, 2 columns

% Selected indices for histogram
selected_indices = ellipse_indices;

% Auto color generation
colors = lines(5);

% Pure component (G, S)
G1 = G_on_long;
G2 = G_on_short;
S1 = S_on_long;
S2 = S_on_short;

% Build weight matrix A (3x2 components mapped from (G,S,1))
A = [G1, G2; S1, S2; 1, 1];

% Define W range
W_range = linspace(0, 1, 500); 

% W1 histogram
nexttile;
hold on;
for idx = 1:length(selected_indices)
    numIndex = selected_indices(idx);
    
    G_values_filtered = G_values_filtered_all{numIndex + 1};
    S_values_filtered = S_values_filtered_all{numIndex + 1};

    W1_filtered = zeros(size(G_values_filtered));
    W2_filtered = zeros(size(G_values_filtered));

    for i = 1:length(G_values_filtered)
        B = [G_values_filtered(i); S_values_filtered(i); 1];
        w = A \ B;

        W1_filtered(i) = min(max(w(1), 0), 1);
        W2_filtered(i) = min(max(w(2), 0), 1);
    end

    mu_W1 = mean(W1_filtered);
    sigma_W1 = std(W1_filtered);
    pdf_W1 = normpdf(W_range, mu_W1, sigma_W1);

    histogram(W1_filtered, 10, 'Normalization', 'pdf', 'FaceAlpha', 0.5, ...
        'EdgeColor', 'none', 'FaceColor', colors(idx, :));

    plot(W_range, pdf_W1, '--', 'Color', colors(idx, :) * 0.8, 'LineWidth', 2);
end

title('Weight Distribution of Pure Component 1 (W_1)', 'FontSize', 12);
xlabel('W_1', 'FontSize', 10);
ylabel('Probability Density', 'FontSize', 10);
grid on;
hold off;

% W2 histogram
nexttile;
hold on;
for idx = 1:length(selected_indices)
    numIndex = selected_indices(idx);
    
    G_values_filtered = G_values_filtered_all{numIndex + 1};
    S_values_filtered = S_values_filtered_all{numIndex + 1};

    W1_filtered = zeros(size(G_values_filtered));
    W2_filtered = zeros(size(G_values_filtered));

    for i = 1:length(G_values_filtered)
        B = [G_values_filtered(i); S_values_filtered(i); 1];
        w = A \ B;

        W1_filtered(i) = min(max(w(1), 0), 1);
        W2_filtered(i) = min(max(w(2), 0), 1);
    end

    mu_W2 = mean(W2_filtered);
    sigma_W2 = std(W2_filtered);
    pdf_W2 = normpdf(W_range, mu_W2, sigma_W2);

    histogram(W2_filtered, 10, 'Normalization', 'pdf', 'FaceAlpha', 0.5, ...
        'EdgeColor', 'none', 'FaceColor', colors(idx, :));

    plot(W_range, pdf_W2, '--', 'Color', colors(idx, :) * 0.8, 'LineWidth', 2);
end

title('Weight Distribution of Pure Component 2 (W_2)', 'FontSize', 12);
xlabel('W_2', 'FontSize', 10);
ylabel('Probability Density', 'FontSize', 10);
grid on;
hold off;
