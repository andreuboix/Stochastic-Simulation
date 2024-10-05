function M = tune_scaling_constant(c, xi)
    % Range for x 
    seed = 1234;
    rng(seed);
    x_values = linspace(1e-3, 1e6, 1e7);  % the pdf is not defined on x = 0, so we avoid it.
    
    % Compute f(x) (Inverse Gaussian PDF)
    f_values = arrayfun(@(x) inverse_gaussian_pdf(x, c, xi), x_values);
    
    % Compute g(x) (Exponential PDF)
    lambda = xi^2/2;
    g_values = arrayfun(@(x) exponential_proposal(x, lambda), x_values);
    
    % Compute the ratio f(x) / g(x) for values above the tolerance
    ratio_values = f_values ./ g_values;
    
    M = max(ratio_values);
    
    fprintf('Tuned scaling constant M: %.4f\n', M);
    
    figure;
    plot(x_values, f_values, 'r', 'LineWidth', 2); % Plot f(x) in red
    hold on;
    plot(x_values, M * g_values, 'b--', 'LineWidth', 2); % Plot M * g(x) in blue
    xlabel('x');
    ylabel('Density');
    title('Comparison of f(x) and C * g(x)');
    legend('f(x) (Inverse Gaussian PDF)', 'C * g(x) (Scaled Exponential PDF)', 'Location', 'Best');
    grid on;
    hold off;
end

% Inverse Gaussian PDF
function pdf = inverse_gaussian_pdf(x, c, xi)
    if x <= 0
        pdf = 0;
    else
        pdf = (c / (x^(3/2) * sqrt(2 * pi))) * exp(xi * c - 0.5 * (c^2 / x + xi^2 * x));
    end
end

% Exponential proposal PDF
function pdf = exponential_proposal(x, lambda)
    if x <= 0
        pdf = 0;
    else
        pdf = lambda * exp(-lambda * x);
    end
end
