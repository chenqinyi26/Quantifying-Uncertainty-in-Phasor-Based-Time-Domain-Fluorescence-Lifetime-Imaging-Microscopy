clc;
clear all;

%% Read multiple CSV files
filepaths = {
   'laser1.csv',
   'laser2.csv',
   'laser5.csv',
};

numFiles = length(filepaths);
G_all = cell(numFiles, 1);
S_all = cell(numFiles, 1);
total_photons_all = cell(numFiles, 1); % Store total photon count per file

for f = 1:numFiles
    % Read CSV file
    opts = detectImportOptions(filepaths{f}, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
    opts.DataLines = [3 inf];  % Start reading from line 3
    data = table2array(readtable(filepaths{f}, opts));

    % Extract raw histogram data
    histAll_raw = data(:, 2:2:800)'; % 40 rows matrix

    % Remove NaN values
    nan_col = find(any(isnan(histAll_raw), 1), 1, 'first');
    if ~isempty(nan_col)
        histAll_raw = histAll_raw(:, 1:nan_col-1);
    end

    % Extract time axis
    time = data(1:size(histAll_raw, 2), 1);
    dt_all = diff(time);
    dt = mean(dt_all) * 1e-9;

    % Estimate background noise
    first_zero_idx = find(data(:, 802) == 0, 1, 'first');
    if ~isempty(first_zero_idx)
        laser_time = data(first_zero_idx, 801);
        time_bgnoise = find(time < laser_time, 1, 'last');
        third_start = round(time_bgnoise / 5);
        third_end = round(2 * time_bgnoise / 5);
        histAll_bgnoise = mean(histAll_raw(:, third_start:third_end), 2);
        histAll_cleaned = histAll_raw - histAll_bgnoise;
    end

    % Find last valid time from column 801
    valid_times = data(~isnan(data(:, 801)) & data(:, 801) ~= 0, 801);
    cutoff_time = valid_times(end);

    % Find index where time exceeds cutoff_time
    time_start = find(time > cutoff_time, 1, 'first');

    % Keep only valid decay portion after cutoff_time
    histAll_cleaned = histAll_cleaned(:, time_start:end);
    time = time(time_start:end);

    % Remove negative values
    histAll_cleaned(histAll_cleaned < 0) = 0;

    % Calculate G, S
    numSelectedBins = 40;
    StartBin = 5;
    t_end = dt * numSelectedBins;
    f_val = 1 / t_end;
    omega = 2 * pi * f_val;
    t = dt / 2 : dt : t_end - (dt / 2);

    if size(histAll_cleaned, 2) > (StartBin + numSelectedBins)
        histAll_selected = histAll_cleaned(:, StartBin+1:StartBin+numSelectedBins);
    else
        histAll_selected = histAll_cleaned;
    end

    numPixels = size(histAll_selected, 1);
    G_values = zeros(numPixels, 1);
    S_values = zeros(numPixels, 1);
    total_photons_all{f} = sum(histAll_selected, 2); % Total photon count per pixel

    for i = 1:numPixels
        total_photons = sum(histAll_selected(i, :));
        G_numerator = sum(histAll_selected(i, :) .* cos(omega * t));
        S_numerator = sum(histAll_selected(i, :) .* sin(omega * t));

        if total_photons > 0
            G_values(i) = G_numerator / total_photons;
            S_values(i) = S_numerator / total_photons;
        else
            G_values(i) = NaN;
            S_values(i) = NaN;
        end
    end

    G_all{f} = G_values;
    S_all{f} = S_values;
end

%% Initialize storage variables
G_means = cell(numFiles, 1);
S_means = cell(numFiles, 1);
Var_G_all = cell(numFiles, 1);
Var_S_all = cell(numFiles, 1);
Cov_GS_all = cell(numFiles, 1);

%% Plot G-S 2D histogram (Phasor plot) with 3σ uncertainty ellipses and reference semicircle
figure;
numCols = ceil(sqrt(numFiles)); % Number of subplot columns
numRows = ceil(numFiles / numCols); % Number of subplot rows
theta = linspace(0, 2*pi, 100); % Angle points for ellipse
n_std = 3; % 3σ error ellipse
N0 = 1e5; % Initial photon number
fixed_xlim = [0.24, 0.3]; % Unified G-axis range
fixed_ylim = [0.42, 0.48]; % Unified S-axis range

% Reference semicircle (only in visible range)
G_ref = linspace(fixed_xlim(1), fixed_xlim(2), 100);
S_ref = sqrt(0.25 - (G_ref - 0.5).^2);
valid_idx = S_ref >= fixed_ylim(1) & S_ref <= fixed_ylim(2);
G_ref = G_ref(valid_idx);
S_ref = S_ref(valid_idx);

for f = 1:numFiles
    subplot(numRows, numCols, f);
    
    % Set fixed axis limits
    xlim(fixed_xlim);
    ylim(fixed_ylim);
    axis equal;

    % 2D histogram
    histogram2(G_all{f}, S_all{f}, [10, 10], 'Normalization', 'count', ...
           'DisplayStyle', 'tile', 'EdgeColor', 'none');
    colorbar;

    % Calculate uncertainty from theoretical Nbox
    G_mean = mean(G_all{f}, 'omitnan');
    S_mean = mean(S_all{f}, 'omitnan');
    mean_GS = [G_mean; S_mean];

    % Store means
    G_means{f} = G_mean;
    S_means{f} = S_mean;

    % Estimate tau from G, S
    tau_mean = S_mean / (G_mean * omega);

    % Simulated mono-exponential decay curve
    component_curve = N0 * exp(-t / tau_mean);
    total_photon_component = sum(component_curve);

    % Scale to match experimental photon count
    total_photons = mean(total_photons_all{f}, 'omitnan');
    normalized_P = total_photons / total_photon_component;
    Nbox = component_curve * normalized_P;

    % Compute variance and covariance
    Var_G = sum(Nbox .* (cos(omega * t) - G_mean).^2, 'omitnan') / (sum(Nbox)^2);
    Var_S = sum(Nbox .* (sin(omega * t) - S_mean).^2, 'omitnan') / (sum(Nbox)^2);
    Cov_GS = sum(Nbox .* (cos(omega * t) - G_mean) .* (sin(omega * t) - S_mean), 'omitnan') / (sum(Nbox)^2);

    % Store results
    Var_G_all{f} = Var_G;
    Var_S_all{f} = Var_S;
    Cov_GS_all{f} = Cov_GS;

    % Form covariance matrix
    cov_GS_matrix = [Var_G, Cov_GS; Cov_GS, Var_S];

    % Eigendecomposition for ellipse
    [eig_vec, eig_val] = eig(cov_GS_matrix);
    radii = n_std * sqrt(diag(eig_val));
    ellipse = [cos(theta); sin(theta)] .* radii;
    rotated_ellipse = eig_vec * ellipse;

    % Translate to mean
    ellipse_x = rotated_ellipse(1, :) + mean_GS(1);
    ellipse_y = rotated_ellipse(2, :) + mean_GS(2);

    % Plot ellipse and reference semicircle
    hold on;
    plot(ellipse_x, ellipse_y, 'r', 'LineWidth', 2);
    plot(G_ref, S_ref, 'black--', 'LineWidth', 2);
    hold off;

    % Fix axis again
    xlim(fixed_xlim);
    ylim(fixed_ylim);
    
    xlabel('G Value');
    ylabel('S Value');
    title(['2D Histogram (Dataset ', num2str(f), ')']);
    grid on;
end

%% Plot histograms of G and S with theoretical uncertainty curves
figure;

numBins = 30; % Number of histogram bins
colors = lines(numFiles); % Different colors for datasets

% Determine bin edges based on min/max of all datasets
G_min = min(cellfun(@min, G_all));
G_max = max(cellfun(@max, G_all));
S_min = min(cellfun(@min, S_all));
S_max = max(cellfun(@max, S_all));

G_edges = linspace(G_min, G_max, numBins + 1);
S_edges = linspace(S_min, S_max, numBins + 1);

% Compute overall means
G_all_concat = vertcat(G_all{:});
G_mean = mean(G_all_concat);

S_all_concat = vertcat(S_all{:});
S_mean = mean(S_all_concat);

% Left: G histogram
subplot(1, 2, 1);
hold on;
for f = 1:numFiles
    [G_counts, ~] = histcounts(G_all{f}, G_edges);
    G_bin_centers = (G_edges(1:end-1) + G_edges(2:end)) / 2;
    bin_width = G_edges(2) - G_edges(1);

    bar(G_bin_centers, G_counts, 'FaceColor', colors(f, :), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceAlpha', 0.5);

    G_std = sqrt(Var_G_all{f});
    if G_std > 0
        fine_x = linspace(G_min, G_max, 500);
        G_theoretical = numel(G_all{f}) * bin_width * (1 / (sqrt(2 * pi) * G_std)) .* exp(-((fine_x - G_mean).^2) ./ (2 * Var_G_all{f}));
        plot(fine_x, G_theoretical, '--', 'Color', colors(f, :), 'LineWidth', 2);
    end
end
xlabel('G Value');
ylabel('Counts');
title('Histogram of G Values with Uncertainty Model');
grid on;
xlim([0.25 0.3]);
ylim([0 80]);
box on;
hold off;

% Right: S histogram
subplot(1, 2, 2);
hold on;
for f = 1:numFiles
    [S_counts, ~] = histcounts(S_all{f}, S_edges);
    S_bin_centers = (S_edges(1:end-1) + S_edges(2:end)) / 2;
    bin_width = S_edges(2) - S_edges(1);

    bar(S_bin_centers, S_counts, 'FaceColor', colors(f, :), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceAlpha', 0.5);

    S_std = sqrt(Var_S_all{f});
    if S_std > 0
        fine_y = linspace(S_min, S_max, 500);
        S_theoretical = numel(S_all{f}) * bin_width * (1 / (sqrt(2 * pi) * S_std)) .* exp(-((fine_y - S_mean).^2) ./ (2 * Var_S_all{f}));
        plot(fine_y, S_theoretical, '--', 'Color', colors(f, :), 'LineWidth', 2);
    end
end
xlabel('S Value');
ylabel('Counts');
title('Histogram of S Values with Uncertainty Model');
grid on;
xlim([0.42 0.47]);
ylim([0 80]);
box on;
hold off;
