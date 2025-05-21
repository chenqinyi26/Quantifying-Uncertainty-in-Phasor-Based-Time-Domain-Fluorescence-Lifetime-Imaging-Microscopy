function [Cov_matrix, G, S] = uncertainty_fluorescent_dyes(total_photons, P1, P2, mean_Tau_purelong, mean_Tau_pureshort)
    % uncertainty_fluorescent_dyes - Calculate the covariance matrix for G and S
    % given component weights (P1, P2) and total photon count.
    %
    % Inputs:
    %   total_photons         - Total number of photons
    %   P1                    - Weight of component 1
    %   P2                    - Weight of component 2
    %   mean_Tau_purelong     - Lifetime of the long component (tau1)
    %   mean_Tau_pureshort    - Lifetime of the short component (tau2)
    %
    % Outputs:
    %   Cov_matrix            - 2x2 covariance matrix of [G, S]
    %   G                     - Computed G value of the mixture
    %   S                     - Computed S value of the mixture

    % Load global variables
    global dt t_end f omega t numSelectedBins

    N0 = 11e1;  % Initial photon count (arbitrary scaling factor)

    % Assign lifetimes
    tau1 = mean_Tau_purelong;
    tau2 = mean_Tau_pureshort;

    % Generate decay curves for each component
    component1_curve = N0 * exp(-t / tau1);
    component2_curve = N0 * exp(-t / tau2);

    % Compute total photons per component
    total_photon_component1 = sum(component1_curve);
    total_photon_component2 = sum(component2_curve);

    % Normalize weights so that the total photon count is preserved
    total_photons_components = total_photon_component1 + total_photon_component2;
    normalized_total = total_photons / total_photons_components;

    normalized_P1 = P1 * total_photons_components / total_photon_component1;
    normalized_P2 = P2 * total_photons_components / total_photon_component2;

    % Scale decay curves to match photon proportions
    adjusted_component1 = component1_curve * normalized_P1 * normalized_total;
    adjusted_component2 = component2_curve * normalized_P2 * normalized_total;

    % Create the final mixed decay curve
    mixed_decay_curve = adjusted_component1 + adjusted_component2;

    % Compute G and S for the mixed curve
    photon_counts = mixed_decay_curve;
    total_photons = sum(photon_counts);

    G_numerator = sum(photon_counts .* cos(omega * t));
    S_numerator = sum(photon_counts .* sin(omega * t));

    G = G_numerator / total_photons;
    S = S_numerator / total_photons;

    % Compute variance and covariance
    N_box = photon_counts;
    cos_values = cos(omega * t);
    sin_values = sin(omega * t);

    numerator_G = sum(N_box .* (cos_values - G).^2);
    numerator_S = sum(N_box .* (sin_values - S).^2);
    denominator = (sum(N_box))^2;

    Var_G_noisy = numerator_G / denominator;
    Var_S_noisy = numerator_S / denominator;

    numerator_cov_GS = sum(N_box .* (cos_values - G) .* (sin_values - S));
    Cov_GS = numerator_cov_GS / denominator;

    % Construct covariance matrix
    Cov_matrix = [Var_G_noisy, Cov_GS;
                  Cov_GS,      Var_S_noisy];
end

