clear;
clc;

%% Initialize parameters
dt = 1e-10;               % Time bin width
t_end = 5e-9;             % Total acquisition time
f = 1/t_end;              % Frequency
n = 1;                    % Harmonic number
omega = 2 * n * pi * f;   % Angular frequency
t = dt/2:dt:t_end-(dt/2); % Time vector
N0 = 1e5;

% Component lifetimes
tau1 = 0.1e-9; 
tau2 = 1e-9;
tau3 = 4e-9;

% Calculate G and S for pure components
G1 = 1 / (1 + (omega * tau1)^2);
S1 = (omega * tau1) / (1 + (omega * tau1)^2);
G2 = 1 / (1 + (omega * tau2)^2);
S2 = (omega * tau2) / (1 + (omega * tau2)^2);
G3 = 1 / (1 + (omega * tau3)^2);
S3 = (omega * tau3) / (1 + (omega * tau3)^2);

% Define different total photon counts
photon_counts_list = [100, 500, 1000, 2000];

% Normalized weights
P1 = 0.2;
P2 = 0.3;
P3 = 0.5;

% Generate decay curves
component1_curve = exp(-t / tau1);
component2_curve = exp(-t / tau2);
component3_curve = exp(-t / tau3);

% Construct coefficient matrix A (3x3)
A = [G1, G2, G3;
     S1, S2, S3;
     1,  1,  1];  

% Color settings
colors = {'b', 'g', 'r', 'c'}; % Colors for each photon count
num_iterations = 10000;       % Monte Carlo simulation times

% Containers to store G, S, weights, std and covariance
G_values = containers.Map();
S_values = containers.Map();
w_values_map = containers.Map();
G_std_map = containers.Map();
S_std_map = containers.Map();
w_std_map = containers.Map();
cov_matrix_map = containers.Map();

for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);

    % Normalize component weights
    total_photon_component1 = sum(component1_curve);
    total_photon_component2 = sum(component2_curve);
    total_photon_component3 = sum(component3_curve);
    
    normalized_P1 = P1 * total_photons / total_photon_component1;
    normalized_P2 = P2 * total_photons / total_photon_component2;
    normalized_P3 = P3 * total_photons / total_photon_component3;
    
    % Compute mixed decay curve
    mixed_decay_curve = normalized_P1 * component1_curve + ...
                        normalized_P2 * component2_curve + ...
                        normalized_P3 * component3_curve;

    % Theoretical G and S
    cos_values = cos(omega * t);
    sin_values = sin(omega * t);
    
    G_theory = sum(mixed_decay_curve .* cos_values) / sum(mixed_decay_curve);
    S_theory = sum(mixed_decay_curve .* sin_values) / sum(mixed_decay_curve);

    % Variance and covariance
    Var_G_noisy = sum(mixed_decay_curve .* (cos_values - G_theory).^2) / (sum(mixed_decay_curve)^2);
    Var_S_noisy = sum(mixed_decay_curve .* (sin_values - S_theory).^2) / (sum(mixed_decay_curve)^2);
    Cov_GS_noisy = sum(mixed_decay_curve .* (cos_values - G_theory) .* (sin_values - S_theory)) / (sum(mixed_decay_curve)^2);

    G_std = sqrt(Var_G_noisy);
    S_std = sqrt(Var_S_noisy);

    G_std_map(num2str(total_photons)) = G_std;
    S_std_map(num2str(total_photons)) = S_std;

    cov_matrix = [Var_G_noisy, Cov_GS_noisy; Cov_GS_noisy, Var_S_noisy];
    cov_matrix_map(num2str(total_photons)) = cov_matrix;

    % Monte Carlo simulations
    G_noisy_values = zeros(1, num_iterations);
    S_noisy_values = zeros(1, num_iterations);
    w_values = zeros(3, num_iterations); 

    for i = 1:num_iterations
        photon_counts_noisy = poissrnd(mixed_decay_curve);
        G_noisy_values(i) = sum(photon_counts_noisy .* cos_values) / sum(photon_counts_noisy);
        S_noisy_values(i) = sum(photon_counts_noisy .* sin_values) / sum(photon_counts_noisy);
        
        % Solve for weights
        b_noisy = [G_noisy_values(i); S_noisy_values(i); 1];
        w_noisy = A \ b_noisy;
        w_values(:, i) = w_noisy;
    end

    G_values(num2str(total_photons)) = G_noisy_values;
    S_values(num2str(total_photons)) = S_noisy_values;
    w_values_map(num2str(total_photons)) = w_values;
    
    w_std = std(w_values, 0, 2);
    w_std_map(num2str(total_photons)) = w_std;
end

%% Phasor plots and distributions in a single figure
figure;
set(gcf, 'Color', 'w'); % White background

% Left column: phasor plots with semicircle and 3σ ellipse
for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);
    G_data = G_values(num2str(total_photons));
    S_data = S_values(num2str(total_photons));
    
    subplot(4, 2, 2*idx-1);
    histogram2(G_data, S_data, [30, 30], 'DisplayStyle', 'tile', 'Normalization', 'pdf', 'EdgeColor', 'none');
    colorbar;
    hold on;
    
    theta = linspace(0, pi, 100);
    x = 0.5 + 0.5 * cos(theta);
    y = 0.5 * sin(theta);
    plot(x, y, 'k--', 'LineWidth', 1.5);
    
    plot(G1, S1, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    text(G1 + 0.02, S1, 'G1,S1', 'FontSize', 12, 'Color', 'k');
    plot(G2, S2, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    text(G2 + 0.02, S2, 'G2,S2', 'FontSize', 12, 'Color', 'k');
    plot(G3, S3, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    text(G3 + 0.02, S3, 'G3,S3', 'FontSize', 12, 'Color', 'k');
    plot([G1, G2, G3, G1], [S1, S2, S3, S1], 'k--', 'LineWidth', 1.5);
    
    Cov_matrix = cov_matrix_map(num2str(total_photons));
    [V, D] = eig(Cov_matrix);
    eigvals = diag(D);
    width = 3 * sqrt(eigvals(1));
    height = 3 * sqrt(eigvals(2));
    angle = atan2(V(2,1), V(1,1)) * (180/pi);
    
    theta = linspace(0, 2*pi, 100);
    ellipse_x = width * cos(theta);
    ellipse_y = height * sin(theta);
    R = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
    ellipse_rotated = R * [ellipse_x; ellipse_y];
    
    plot(mean(G_data) + ellipse_rotated(1, :), mean(S_data) + ellipse_rotated(2, :), 'r', 'LineWidth', 2);
    
    xlabel('G', 'FontSize', 12, 'Color', 'k');
    ylabel('S', 'FontSize', 12, 'Color', 'k');
    title(sprintf('%d Photons Phasor Plot', total_photons), 'FontSize', 12, 'Color', 'k');
    set(gca, 'XColor', 'k', 'YColor', 'k');
    axis equal;
    grid on;
end

% Distribution of G and S
distributions = {'G', 'S', 'w1'};
data_maps = {G_values, S_values};
std_maps = {G_std_map, S_std_map};

w1_data = w_values(1, :);
w2_data = w_values(2, :);
w3_data = w_values(3, :);
theoretical_values = {G_theory, S_theory, mean(w1_data)};
num_bins = 100;
colors = {'k', 'b', 'g', 'r'};
line_styles = {'--', '-.', ':'};

for i = 1:2
    all_data = [];
    for idx = 1:length(photon_counts_list)
        total_photons = photon_counts_list(idx);
        data = data_maps{i}(num2str(total_photons));
        all_data = [all_data; data(:)];
    end
    bin_edges = linspace(min(all_data), max(all_data), num_bins);
    bin_width = diff(bin_edges(1:2));

    subplot(4, 2, 2 * i);
    hold on;

    for idx = 1:length(photon_counts_list)
        total_photons = photon_counts_list(idx);
        data = data_maps{i}(num2str(total_photons));
        data_std = std_maps{i}(num2str(total_photons));
        histogram(data(:), 'BinEdges', bin_edges, 'FaceColor', colors{idx}, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.5);

        x = linspace(min(all_data), max(all_data), 100);
        y_scale_factor = length(data(:)) * bin_width;
        plot(x, normpdf(x, theoretical_values{i}, data_std) * y_scale_factor, ...
             '--', 'Color', colors{idx}, 'LineWidth', 2);
    end

    title(sprintf('%s Noisy Distribution', distributions{i}), 'FontSize', 12, 'Color', 'k');
    xlim([0.1 0.5]);
    xlabel(distributions{i}, 'FontSize', 12, 'Color', 'k');
    ylabel('Count', 'FontSize', 12, 'Color', 'k');
    set(gca, 'XColor', 'k', 'YColor', 'k');
    box on;
    grid on;
end

% Distribution of w1, w2, w3
subplot(4, 2, 6);
hold on;

all_w_data = [];
for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);
    w_values = w_values_map(num2str(total_photons));
    all_w_data = [all_w_data, w_values(:)];
end
all_w_data = all_w_data(:);
bin_edges = linspace(min(all_w_data), max(all_w_data), 2*num_bins);
bin_width = diff(bin_edges(1:2));

for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);
    w_values = w_values_map(num2str(total_photons));

    for j = 1:3
        data = w_values(j, :);
        color = colors{idx};
        histogram(data, 'BinEdges', bin_edges, 'EdgeColor', 'none', ...
                  'FaceColor', color, 'FaceAlpha', 0.3);

        x = linspace(min(all_w_data), max(all_w_data), 100);
        y_scale_factor = length(data) * bin_width;
        plot(x, normpdf(x, mean(data), std(data)) * y_scale_factor, ...
             line_styles{j}, 'Color', color, 'LineWidth', 2);
    end
end

xlim([0 0.8]);
title('w1, w2, w3 Noisy Distribution', 'FontSize', 12, 'Color', 'k');
xlabel('Weight Value', 'FontSize', 12, 'Color', 'k');
ylabel('Count', 'FontSize', 12, 'Color', 'k');
box on;
grid on;

% Apply global text and axis settings
allText = findall(gcf, 'Type', 'text');
for i = 1:length(allText)
    set(allText(i), 'Color', 'k');
end
allAxes = findall(gcf, 'Type', 'axes');
for i = 1:length(allAxes)
    set(allAxes(i), 'XColor', 'k', 'YColor', 'k', 'FontName', 'Arial');
end
hold off;

%% Monte Carlo point percentage within 1σ, 2σ, 3σ ellipses
sigma_levels = [1, 2, 3];
percentages = zeros(length(photon_counts_list), length(sigma_levels));

for idx = 1:length(photon_counts_list)
    total_photons = photon_counts_list(idx);
    G_sim = G_values(num2str(total_photons));
    S_sim = S_values(num2str(total_photons));
    G_mean = mean(G_sim);
    S_mean = mean(S_sim);
    Cov_matrix = cov_matrix_map(num2str(total_photons));
    Cov_inv = inv(Cov_matrix);
    
    num_points = length(G_sim);
    distances = zeros(1, num_points);
    
    for i = 1:num_points
        delta_x = [G_sim(i) - G_mean; S_sim(i) - S_mean];
        distances(i) = delta_x' * Cov_inv * delta_x;
    end
    
    for j = 1:length(sigma_levels)
        sigma_threshold = sigma_levels(j)^2;
        percentages(idx, j) = sum(distances <= sigma_threshold) / num_points * 100;
    end
end

% Display results
disp('Percentage of Monte Carlo points within 1σ, 2σ, 3σ:');
fprintf('Photons\t 1σ(%%)\t 2σ(%%)\t 3σ(%%)\n');
for idx = 1:length(photon_counts_list)
    fprintf('%d\t %.2f\t %.2f\t %.2f\n', photon_counts_list(idx), ...
            percentages(idx, 1), percentages(idx, 2), percentages(idx, 3));
end

% Subplot 8: plot percentages vs. photon count
subplot(4, 2, 8);
hold on;
sigma_colors = [1 0 0; 0 1 0; 0 0 1]; % RGB for 1σ, 2σ, 3σ
markers = {'o', 's', 'd'};

for j = 1:length(sigma_levels)
    plot(photon_counts_list, percentages(:, j), ...
         '-o', 'LineWidth', 2, ...
         'Marker', markers{j}, ...
         'Color', sigma_colors(j, :), ...
         'MarkerFaceColor', sigma_colors(j, :), ...
         'DisplayName', sprintf('%dσ', sigma_levels(j)));
end

xlabel('Photon Count', 'FontSize', 12, 'Color', 'k');
ylabel('Percentage within σ (%)', 'FontSize', 12, 'Color', 'k');
ylim([30 110]);
title('Percentage of Points within 1σ, 2σ, 3σ', 'FontSize', 12, 'Color', 'k');
legend('Location', 'southeast');
set(gca, 'XColor', 'k', 'YColor', 'k', 'FontName', 'Arial');
grid on;
box on;
