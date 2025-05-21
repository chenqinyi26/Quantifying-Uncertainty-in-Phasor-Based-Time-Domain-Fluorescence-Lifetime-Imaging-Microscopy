clear;
clc;

%% Initialize parameters
dt = 1e-10;
t_end = 5e-9;
f = 1 / t_end;
n = 1;
omega = 2 * n * pi * f;
t = dt/2:dt:t_end-(dt/2);
N0 = 1e5;

% Decay lifetimes (two components)
tau1 = 0.5e-9;
tau2 = 4e-9;

% Phasor coordinates (G, S)
G1 = 1 / (1 + (omega * tau1)^2);
S1 = (omega * tau1) / (1 + (omega * tau1)^2);
G2 = 1 / (1 + (omega * tau2)^2);
S2 = (omega * tau2) / (1 + (omega * tau2)^2);

% Photon counts to simulate
photon_counts_list = [100, 500, 1000, 2000];

% Weight proportions (two components)
P1 = 0.5;
P2 = 0.5;

% Decay curves
component1_curve = exp(-t / tau1);
component2_curve = exp(-t / tau2);

% Coefficient matrix A (2 components)
A = [G1, G2;
     S1, S2;
     1,  1];

colors = {'b', 'g', 'r', 'c'};
num_iterations = 10000;

% Data storage maps
G_values = containers.Map();
S_values = containers.Map();
w_values_map = containers.Map();
G_std_map = containers.Map();
S_std_map = containers.Map();
w_std_map = containers.Map();
cov_matrix_map = containers.Map();

% Monte Carlo simulation for each photon count
for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);
    
    % Normalize weights
    total1 = sum(component1_curve);
    total2 = sum(component2_curve);
    norm_P1 = P1 * total_photons / total1;
    norm_P2 = P2 * total_photons / total2;

    % Mixed decay curve
    mixed_curve = norm_P1 * component1_curve + ...
                  norm_P2 * component2_curve;

    % Theoretical G and S
    cos_vals = cos(omega * t);
    sin_vals = sin(omega * t);
    G_theory = sum(mixed_curve .* cos_vals) / sum(mixed_curve);
    S_theory = sum(mixed_curve .* sin_vals) / sum(mixed_curve);

    % Variance and covariance (uncertainty model)
    Var_G = sum(mixed_curve .* (cos_vals - G_theory).^2) / (sum(mixed_curve)^2);
    Var_S = sum(mixed_curve .* (sin_vals - S_theory).^2) / (sum(mixed_curve)^2);
    Cov_GS = sum(mixed_curve .* (cos_vals - G_theory) .* (sin_vals - S_theory)) / (sum(mixed_curve)^2);

    % Store uncertainty values
    G_std_map(num2str(total_photons)) = sqrt(Var_G);
    S_std_map(num2str(total_photons)) = sqrt(Var_S);
    cov_matrix_map(num2str(total_photons)) = [Var_G, Cov_GS; Cov_GS, Var_S];

    % Monte Carlo sampling
    G_noisy = zeros(1, num_iterations);
    S_noisy = zeros(1, num_iterations);
    w_values = zeros(2, num_iterations);

    for i = 1:num_iterations
        photons = poissrnd(mixed_curve);
        G_noisy(i) = sum(photons .* cos_vals) / sum(photons);
        S_noisy(i) = sum(photons .* sin_vals) / sum(photons);
        b = [G_noisy(i); S_noisy(i); 1];
        w_values(:, i) = A \ b;
    end

    % Store results
    G_values(num2str(total_photons)) = G_noisy;
    S_values(num2str(total_photons)) = S_noisy;
    w_values_map(num2str(total_photons)) = w_values;
    w_std_map(num2str(total_photons)) = std(w_values, 0, 2);
end

%% Plot phasor and histogram panels: left = phasor, right = distributions
figure;
set(gcf, 'Color', 'w');

% Phasor plots
for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);
    G_data = G_values(num2str(total_photons));
    S_data = S_values(num2str(total_photons));

    subplot(4, 2, 2*idx-1);
    histogram2(G_data, S_data, [30, 30], 'DisplayStyle', 'tile', ...
               'Normalization', 'pdf', 'EdgeColor', 'none');
    colorbar;
    hold on;

    % Semicircle
    theta = linspace(0, pi, 100);
    x = 0.5 + 0.5 * cos(theta);
    y = 0.5 * sin(theta);
    plot(x, y, 'k--', 'LineWidth', 1.5);

    % Component phasors and connection
    plot(G1, S1, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    text(G1 + 0.02, S1, 'G1,S1', 'FontSize', 12, 'Color', 'k');
    plot(G2, S2, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    text(G2 + 0.02, S2, 'G2,S2', 'FontSize', 12, 'Color', 'k');
    plot([G1, G2], [S1, S2], 'k--', 'LineWidth', 1.5);

    % 3σ ellipse
    Cov_matrix = cov_matrix_map(num2str(total_photons));
    [V, D] = eig(Cov_matrix);
    width = 3 * sqrt(D(1,1));
    height = 3 * sqrt(D(2,2));
    angle = atan2(V(2,1), V(1,1)) * 180/pi;
    theta = linspace(0, 2*pi, 100);
    ellipse = [width * cos(theta); height * sin(theta)];
    R = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
    ellipse = R * ellipse;
    plot(mean(G_data) + ellipse(1,:), mean(S_data) + ellipse(2,:), 'r', 'LineWidth', 2);

    xlabel('G'); ylabel('S');
    title(sprintf('%d Photons Phasor Plot', total_photons));
    axis equal; grid on;
end

%% G and S histograms + normal fits
distributions = {'G', 'S'};
data_maps = {G_values, S_values};
std_maps = {G_std_map, S_std_map};
theoretical_values = {G_theory, S_theory};
colors = {'k', 'b', 'g', 'r'};
num_bins = 100;

for i = 1:2
    all_data = [];
    for idx = 1:length(photon_counts_list)
        data = data_maps{i}(num2str(photon_counts_list(idx)));
        all_data = [all_data; data(:)];
    end
    bin_edges = linspace(min(all_data), max(all_data), num_bins + 5);
    bin_width = diff(bin_edges(1:2));

    subplot(4, 2, 2 * i);
    hold on;
    for idx = 1:length(photon_counts_list)
        photons = photon_counts_list(idx);
        data = data_maps{i}(num2str(photons));
        data_std = std_maps{i}(num2str(photons));
        histogram(data(:), 'BinEdges', bin_edges, 'FaceColor', colors{idx}, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.5);
        x = linspace(min(all_data), max(all_data), 100);
        y = normpdf(x, theoretical_values{i}, data_std) * length(data(:)) * bin_width;
        plot(x, y, '--', 'Color', colors{idx}, 'LineWidth', 2);
    end
    title(sprintf('%s Noisy Distribution', distributions{i}));
    xlim([0.15 0.55]);
    xlabel(distributions{i}); ylabel('Count');
    box on; grid on;
end

%% w1 and w2 histograms
subplot(4, 2, 6);
hold on;
all_w_data = [];
for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);
    w_vals = w_values_map(num2str(total_photons));
    all_w_data = [all_w_data, w_vals(1, :), w_vals(2, :)];
end

if isempty(all_w_data)
    error('all_w_data is empty. Run the simulation code first.');
end

bin_edges = linspace(min(all_w_data), max(all_w_data), 2*num_bins);
bin_width = diff(bin_edges(1:2));
line_styles = {'--', '-.'};

for idx = 1:length(photon_counts_list)
    photons = photon_counts_list(idx);
    w = w_values_map(num2str(photons));
    for j = 1:2
        data = w(j, :);
        color = colors{idx};
        histogram(data, 'BinEdges', bin_edges, 'EdgeColor', 'none', ...
                  'FaceColor', color, 'FaceAlpha', 0.3);
        x = linspace(min(all_w_data), max(all_w_data), 100);
        y = normpdf(x, mean(data), std(data)) * length(data) * bin_width;
        plot(x, y, line_styles{j}, 'Color', color, 'LineWidth', 2);
    end
end
title('w1, w2 Noisy Distribution');
xlabel('Weight Value'); ylabel('Count');
xlim([0.2 0.8]);
box on; grid on;

%% σ coverage percentage
sigma_levels = [1, 2, 3];
percentages = zeros(length(photon_counts_list), length(sigma_levels));
for idx = 1:length(photon_counts_list)
    photons = photon_counts_list(idx);
    G_sim = G_values(num2str(photons));
    S_sim = S_values(num2str(photons));
    G_mean = mean(G_sim);
    S_mean = mean(S_sim);
    Cov = cov_matrix_map(num2str(photons));
    Cov_inv = inv(Cov);
    distances = zeros(1, length(G_sim));
    for i = 1:length(G_sim)
        delta = [G_sim(i) - G_mean; S_sim(i) - S_mean];
        distances(i) = delta' * Cov_inv * delta;
    end
    for j = 1:length(sigma_levels)
        threshold = sigma_levels(j)^2;
        percentages(idx, j) = sum(distances <= threshold) / length(G_sim) * 100;
    end
end

disp('Percentage of Monte Carlo points within 1σ, 2σ, 3σ ellipse:');
fprintf('Photons\t 1σ(%%)\t 2σ(%%)\t 3σ(%%)\n');
for idx = 1:length(photon_counts_list)
    fprintf('%d\t %.2f\t %.2f\t %.2f\n', photon_counts_list(idx), ...
            percentages(idx, 1), percentages(idx, 2), percentages(idx, 3));
end

subplot(4, 2, 8);
hold on;
sigma_colors = [1 0 0; 0 1 0; 0 0 1];
markers = {'o', 's', 'd'};
for j = 1:3
    plot(photon_counts_list, percentages(:, j), ...
         '-o', 'LineWidth', 2, ...
         'Marker', markers{j}, ...
         'Color', sigma_colors(j, :), ...
         'MarkerFaceColor', sigma_colors(j, :));
end
xlabel('Photon Count');
ylabel('Percentage within σ (%)');
ylim([30 110]);
title('Percentage of Points within 1σ, 2σ, 3σ');
legend({'1σ', '2σ', '3σ'}, 'Location', 'southeast');
grid on; box on;
