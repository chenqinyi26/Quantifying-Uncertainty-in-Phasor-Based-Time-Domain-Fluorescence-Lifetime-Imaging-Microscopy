clc;
clear;

%% Global parameter setup
global dt t_end f omega t numSelectedBins StartBin

% Time parameters
dt = 9.77645305514160e-11;      % Time bin width
numSelectedBins = 40;          % Number of time bins to select
StartBin = 5;                  % Start after the rising edge (after peak)
t_end = dt * numSelectedBins; % Total acquisition time
f = 1 / t_end;                 % Frequency (40 MHz)
omega = 2 * pi * f;           % Angular frequency
t = dt/2 : dt : t_end - (dt/2); % Time array

% Pure components
G1 = 0.070; S1 = 0.255; tau1 = 2.269e-09;
G2 = 0.051; S2 = 0.220; tau2 = 2.687e-09;

% Define file path
filePath = 'frame20\combine\histAll_combined.mat';
% filePath = 'frame5\seperate\histAll_combined.mat';

% Load .mat file
data = load(filePath);
fieldNames = fieldnames(data);        % Get variable names
histAll = data.(fieldNames{1});       % Extract the main matrix

%% Calculate G and S values
[numRows, numCols] = size(histAll);   % Get matrix dimensions

% Initialize arrays
G_values = [];
S_values = [];
total_photons_values = [];
direct_lifetime_values = [];
skipped_rows = 0; % Count of skipped rows due to indexing errors

for i = 1:numRows
    % Skip rows if index exceeds valid range
    if StartBin + numSelectedBins - 1 > numCols
        skipped_rows = skipped_rows + 1;
        continue;
    end
    
    % Extract intensity for valid time bins
    validRange = StartBin : (StartBin + numSelectedBins - 1);
    I_t = histAll(i, validRange);

    % Calculate G, S, total photons, and direct lifetime
    G = sum(I_t .* cos(omega * t)) / sum(I_t);
    S = sum(I_t .* sin(omega * t)) / sum(I_t);
    total_photons = sum(I_t);
    direct_lifetime = S / (G * omega);

    % Store results
    G_values = [G_values; G];
    S_values = [S_values; S];
    total_photons_values = [total_photons_values; total_photons];
    direct_lifetime_values = [direct_lifetime_values; direct_lifetime];
end

% Plot 2D histogram (phasor plot)
figure;
histogram2(G_values, S_values, [20, 20], 'DisplayStyle', 'tile', 'Normalization', 'pdf');
colorbar;
hold on;

xlabel('G');
ylabel('S');
title('Phasor Plot (2D Histogram)');

% Plot theoretical semicircle (for single-exponential decay)
theta = linspace(0, pi, 100);
x = 0.5 + 0.5 * cos(theta);
y = 0.5 * sin(theta);
plot(x, y, 'k--', 'LineWidth', 1.5);

% Plot pure components
plot(G1, S1, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); % Component 1
plot(G2, S2, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); % Component 2

axis equal;
grid on;
hold off;

%% Compute weights w1 and w2
w1_values = zeros(length(G_values), 1);
w2_values = zeros(length(G_values), 1);

% Construct coefficient matrix A
A = [G1, G2; S1, S2; 1, 1];
w = zeros(2, length(G_values));

for i = 1:length(G_values)
    B = [G_values(i); S_values(i); 1];
    w = A \ B; % Solve linear system
    w1_values(i) = min(max(w(1), 0), 1); % Clamp between [0, 1]
    w2_values(i) = min(max(w(2), 0), 1);
end

%% Compute lifetime for each pixel
lifetime_values = zeros(size(w1_values));
lifetime_values = compute_lifetime(w1_values, w2_values, total_photons_values, tau1, tau2);

%% Reconstruct original 50×50 intensity image
imageSize = 50;
reconstructedImage = zeros(imageSize, imageSize);

for row = 1:imageSize
    startIdx = (row - 1) * imageSize + 1;
    endIdx = row * imageSize;
    reconstructedImage(row, :) = total_photons_values(startIdx:endIdx);
end

% Normalize for consistent display scaling
minVal = min(reconstructedImage(:));
maxVal = max(reconstructedImage(:));

% Show grayscale image
figure;
imagesc(reconstructedImage, [minVal maxVal]);
colormap(gray);
colorbar;
caxis([minVal maxVal]);
axis image;
title('Intensity Image');

%% Generate 50×50 lifetime map
reconstructedImage = zeros(imageSize, imageSize);
lifetimeImage = zeros(imageSize, imageSize);

for row = 1:imageSize
    startIdx = (row - 1) * imageSize + 1;
    endIdx = row * imageSize;
    reconstructedImage(row, :) = total_photons_values(startIdx:endIdx);
    lifetimeImage(row, :) = lifetime_values(startIdx:endIdx);
    % lifetimeImage(row, :) = direct_lifetime_values(startIdx:endIdx); % (optional)
end

% Use parula colormap
minTau = 2.1e-9;
maxTau = 2.8e-9;

% Map lifetime values to colormap
coloredLifetimeImage = ind2rgb(gray2ind(mat2gray(lifetimeImage, [minTau maxTau]), 256), parula(256));

% Normalize intensity to adjust brightness
intensityFactor = reconstructedImage - min(reconstructedImage(:));
intensityFactor = intensityFactor / max(intensityFactor(:));
intensityFactor = intensityFactor .^ 0.5; % Gamma correction

% Combine brightness-adjusted intensity with lifetime color
finalImage = coloredLifetimeImage .* repmat(intensityFactor, [1, 1, 3]);

% Display final overlaid lifetime-intensity image
figure;
imagesc(finalImage);
axis image;
xlabel('X Pixel');
ylabel('Y Pixel');
title('Lifetime-Intensity Overlay Image');

% Set parula colormap with appropriate colorbar scaling
colormap(gca, parula);
colorbar;
caxis([minTau maxTau]);

%% Compute uncertainty of weights for selected pixels
uncertainty_values = compute_uncertainty(G_values, S_values, w1_values, w2_values, total_photons_values, tau1, tau2);
