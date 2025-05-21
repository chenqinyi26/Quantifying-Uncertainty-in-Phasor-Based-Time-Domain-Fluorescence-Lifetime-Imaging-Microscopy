function plot_gaussian_contour(mean_GS, cov_GS, n_std, lineStyle)
    % Eigenvalue decomposition of the covariance matrix
    [eig_vec, eig_val] = eig(cov_GS); % Eigenvectors and eigenvalues
    radii = n_std * sqrt(diag(eig_val)); % Axes lengths for n standard deviations

    % Generate angle points for the ellipse
    theta = linspace(0, 2 * pi, 100);
    ellipse = [cos(theta); sin(theta)] .* radii; % Scale unit circle to an ellipse

    % Rotate the ellipse according to the covariance orientation
    rotated_ellipse = eig_vec * ellipse;

    % Translate the ellipse to the mean location
    ellipse_x = rotated_ellipse(1, :) + mean_GS(1);
    ellipse_y = rotated_ellipse(2, :) + mean_GS(2);

    % Plot the ellipse
    plot(ellipse_x, ellipse_y, lineStyle, 'LineWidth', 1.5);
end
