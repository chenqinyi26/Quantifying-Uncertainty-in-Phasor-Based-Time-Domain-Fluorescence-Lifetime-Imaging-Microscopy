function [uncertainty_values] = compute_uncertainty(G_values, S_values, w1_values, w2_values, total_photons_values, tau1, tau2)

global dt t_end f omega t numSelectedBins

    N0 = 11e1; % Initial photon count (arbitrary unit)

    % Generate decay curves for component 1 and component 2
    component1_curve = N0 * exp(-t / tau1);
    component2_curve = N0 * exp(-t / tau2);

    % Compute total photons for each component (before mixing)
    total_photon_component1 = sum(component1_curve);
    total_photon_component2 = sum(component2_curve);

    % Initialize uncertainty output array
    uncertainty_values = zeros(size(w1_values));

    for i = 1:length(w1_values)
        P1 = w1_values(i);
        P2 = w2_values(i);
        total_photons = total_photons_values(i);

        % Create 100 sample points between 0 and 1 for P1
        if P1 >= 0 && P1 <= 1
            P1_values = linspace(0, 1, 100);
        end
        P2_values = 1 - P1_values;

        probabilities = zeros(1, 100);
        for j = 1:100
            P1_temp = P1_values(j);
            P2_temp = P2_values(j);

            normalized_total = total_photons / (total_photon_component1 + total_photon_component2);
            normalized_P1 = P1_temp * (total_photon_component1 + total_photon_component2) / total_photon_component1;
            normalized_P2 = P2_temp * (total_photon_component1 + total_photon_component2) / total_photon_component2;

            adjusted_component1 = component1_curve * normalized_P1 * normalized_total;
            adjusted_component2 = component2_curve * normalized_P2 * normalized_total;

            mixed_decay_curve = adjusted_component1 + adjusted_component2;

            photon_counts = mixed_decay_curve;
            total_photons = sum(photon_counts);

            G_numerator = sum(photon_counts .* cos(omega * t));
            S_numerator = sum(photon_counts .* sin(omega * t));

            G = G_numerator / total_photons;
            S = S_numerator / total_photons;

            Var_G = sum(photon_counts .* (cos(omega * t) - G).^2) / (sum(photon_counts)^2);
            Var_S = sum(photon_counts .* (sin(omega * t) - S).^2) / (sum(photon_counts)^2);
            Cov_GS = sum(photon_counts .* (cos(omega * t) - G) .* (sin(omega * t) - S)) / (sum(photon_counts)^2);

            C = [Var_G, Cov_GS; Cov_GS, Var_S];

            x = [G_values(i); S_values(i)];
            mu = [G; S];
            p_i = (1 / (2 * pi * sqrt(det(C)))) * exp(-0.5 * (x - mu)' * inv(C) * (x - mu));
            probabilities(j) = p_i;
        end

        % Normalize probabilities
        probabilities = probabilities / sum(probabilities);

        % Compute uncertainty as weighted standard deviation
        uncertainty_values(i) = sqrt(sum(probabilities .* (P1_values - P1).^2));

        % --- Visualization for specific pixel ---
        if i == 2277  % Change to any pixel index of interest
            figure;
            hold on;

            % Normalize probability values for color blending
            probabilities_norm = 0.5 + 0.5 .* (probabilities - min(probabilities)) / (max(probabilities) - min(probabilities));

            % Create low-saturation color maps for P1 (blue) and P2 (red)
            colors_P1 = [ones(length(P1_values), 1) .* (1 - probabilities_norm' * 0.8), ...
                         ones(length(P1_values), 1) .* (1 - probabilities_norm' * 0.8), ...
                         ones(length(P1_values), 1) * 0.8];
            colors_P2 = [ones(length(P2_values), 1) * 0.8, ...
                         ones(length(P2_values), 1) .* (1 - probabilities_norm' * 0.8), ...
                         ones(length(P2_values), 1) .* (1 - probabilities_norm' * 0.8)];

            % Plot weighted probability curves for P1 and P2
            for j = 1:length(P1_values)-1
                plot(P1_values(j:j+1), probabilities(j:j+1) * 100, 'Color', colors_P1(j, :), 'LineWidth', 4);
                plot(P2_values(j:j+1), probabilities(j:j+1) * 100, 'Color', colors_P2(j, :), 'LineWidth', 4);
            end

            % Compute standard deviation from the distribution
            std_dev = sqrt(sum(probabilities .* (P1_values - mean(P1_values)).^2));

            disp(['Pixel ', num2str(i), ' - Standard Deviation: ', num2str(std_dev)]);
            disp(['Pixel ', num2str(i), ' - Total Photons: ', num2str(total_photons)]);

            xlabel('P1 and P2 Values');
            ylabel('Probability Density (%)');
            ylim([0 7]);
            title(['Joint Probability Density of P1 and P2 for Pixel ', num2str(i)]);
            grid on;
            box on;
            hold off;
        end

        % --- Lifetime uncertainty visualization for specific pixel ---
        if i == 2277
            figure;
            hold on;

            % Compute theoretical G and S for pure lifetimes
            G_tau1 = 1 / (1 + (omega * tau1)^2);
            S_tau1 = (omega * tau1) / (1 + (omega * tau1)^2);

            G_tau2 = 1 / (1 + (omega * tau2)^2);
            S_tau2 = (omega * tau2) / (1 + (omega * tau2)^2);

            % Interpolate G, S, and lifetime across 100 points
            G_values_interp = P1_values * G_tau1 + (1 - P1_values) * G_tau2;
            S_values_interp = P1_values * S_tau1 + (1 - P1_values) * S_tau2;
            tau_values = S_values_interp ./ (G_values_interp * omega);

            % Plot fluorescence lifetime vs. probability
            plot(tau_values, probabilities * 100, 'b-', 'LineWidth', 2);

            % Compute mean and standard deviation of lifetime
            tau_mean = sum(probabilities .* tau_values);
            tau_std_dev = sqrt(sum(probabilities .* (tau_values - tau_mean).^2));

            disp(['Pixel ', num2str(i), ' - Fluorescence Lifetime Std Dev: ', num2str(tau_std_dev)]);
            disp(['Pixel ', num2str(i), ' - Total Photons: ', num2str(total_photons)]);

            xlabel('Fluorescence Lifetime (\tau)');
            ylabel('Probability Density (%)');
            title(['Uncertainty Curve of Fluorescence Lifetime for Pixel ', num2str(i)]);
            ylim([0 7]);
            grid on;
            box on;
            hold off;
        end

    end
end

