function [lifetime_values] = compute_lifetime(w1_values, w2_values, total_photons_values, tau1, tau2)
    
    global dt t_end f omega t numSelectedBins
    
    N0 = 11e1; % Initial photon count (arbitrary units)
    
    % Generate decay curves for component 1 and component 2
    component1_curve = N0 * exp(-t / tau1);
    component2_curve = N0 * exp(-t / tau2);
    
    % Compute the total photon count of each component (pre-normalization)
    total_photon_component1 = sum(component1_curve);
    total_photon_component2 = sum(component2_curve);
    
    % Initialize output array for lifetime values
    lifetime_values = zeros(size(w1_values));
    
    for i = 1:length(w1_values)
        P1 = w1_values(i);
        P2 = w2_values(i);
        total_photons = total_photons_values(i);
        
        % Normalize photon contribution for mixture synthesis
        total_photons_components = total_photon_component1 + total_photon_component2;
        normalized_total = total_photons / total_photons_components;
        normalized_P1 = P1 * total_photons_components / total_photon_component1;
        normalized_P2 = P2 * total_photons_components / total_photon_component2;
        
        % Scale the decay curves by their normalized weights
        adjusted_component1 = component1_curve * normalized_P1 * normalized_total;
        adjusted_component2 = component2_curve * normalized_P2 * normalized_total;
        
        % Create the mixed decay curve
        mixed_decay_curve = adjusted_component1 + adjusted_component2;
        
        % Compute G and S phasor coordinates
        photon_counts = mixed_decay_curve;
        total_photons = sum(photon_counts);
        
        G_numerator = sum(photon_counts .* cos(omega * t));
        S_numerator = sum(photon_counts .* sin(omega * t));
        
        G = G_numerator / total_photons;
        S = S_numerator / total_photons;
        
        % Compute fluorescence lifetime from G and S
        lifetime_values(i) = S / (G * omega);
    end
end
