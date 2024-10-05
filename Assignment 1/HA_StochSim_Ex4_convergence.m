function convergence_Ex4
    seed = 1234;
    rng(seed);
    % Parameters of the Inverse Gaussian distribution
    c = 1; 
    xi = 1;
    tolerance = 1e-4; % Tolerance for comparison

    % Initialize parameters for sample generation
    lambda = xi^2/2;  % Parameter for the exponential distribution
    C = 2.4944; % Tuned scaling constant

    % Initialize sample array
    num_samples = 256e6; % Starting number of samples
    samples = zeros(num_samples, 1);
    sample_mean = 0;
    sample_squared_mean = 0;
    sample_variance = 0;
    diff_mean = inf;
    diff_squared_mean = inf;
    diff_variance = inf;

    % Compare with theoretical values
    theoretical_mean = c / xi;
    theoretical_variance = c / xi^3;
    theoretical_squared_mean = theoretical_mean^2 + theoretical_variance;
    
    while (diff_mean > tolerance || diff_squared_mean > tolerance || diff_variance > tolerance) && num_samples < 4e8 
        num_accepted = 0;

        while num_accepted < num_samples
            U1 = rand;
            X_star = -log(U1) / lambda; % inversion method
            U = rand;
            accept_prob = inverse_gaussian_pdf(X_star, c, xi) / (C * exponential_proposal(X_star, lambda));
            if U <= accept_prob
                num_accepted = num_accepted + 1;
                samples(num_accepted) = X_star;
            end
        end
        
        % Compute the sample mean, squared mean, and variance
        sample_mean = mean(samples);
        sample_squared_mean = mean(samples.^2);
        sample_variance = var(samples);
        
        % Compute differences from theoretical values
        diff_mean = abs(sample_mean - theoretical_mean)
        diff_squared_mean = abs(sample_squared_mean - theoretical_squared_mean)
        diff_variance = abs(sample_variance - theoretical_variance)

        % Compute 95% confidence intervals
        z_975 = 1.96; % Z-value for 95% confidence interval
        n = num_samples;
    
        % Confidence interval for mean
        std_error_mean = sqrt(sample_variance / n);
        CI_mean = [sample_mean - z_975 * std_error_mean, sample_mean + z_975 * std_error_mean];
    
        % Confidence interval for squared mean
        sample_variance_X2 = var(samples.^2);
        std_error_squared_mean = sqrt(sample_variance_X2 / n);
        CI_squared_mean = [sample_squared_mean - z_975 * std_error_squared_mean, sample_squared_mean + z_975 * std_error_squared_mean];
    
        % Confidence interval for variance
        chi2_975 = chi2inv(0.975, n-1);
        chi2_025 = chi2inv(0.025, n-1);
        CI_variance = [(n-1) * sample_variance / chi2_975, (n-1) * sample_variance / chi2_025];
    
        % Print the results
        fprintf('Sample Mean: %.4f\n', sample_mean);
        fprintf('95%% Confidence Interval for Mean: [%.4f, %.4f]\n', CI_mean(1), CI_mean(2));
        
        fprintf('Sample Squared Mean: %.4f\n', sample_squared_mean);
        fprintf('95%% Confidence Interval for Squared Mean: [%.4f, %.4f]\n', CI_squared_mean(1), CI_squared_mean(2));
        
        fprintf('Sample Variance: %.4f\n', sample_variance);
        fprintf('95%% Confidence Interval for Variance: [%.4f, %.4f]\n', CI_variance(1), CI_variance(2));

        fprintf('Theoretical Mean: %.4f\n', theoretical_mean);
        fprintf('Theoretical Squared Mean: %.4f\n', theoretical_squared_mean);
        fprintf('Theoretical Variance: %.4f\n', theoretical_variance);
        % Increase number of samples
        if diff_mean > tolerance || diff_squared_mean > tolerance || diff_variance > tolerance
            num_samples = num_samples * 2;
            samples = zeros(num_samples, 1); % Reset samples array
            fprintf('Increasing number of samples to %d\n', num_samples);
        end
        if num_samples > 4e8 && (diff_mean > tolerance || diff_squared_mean > tolerance || diff_variance > tolerance)
            fprintf('Reached the limit for the number of samples\n');
        end
    end
end

% Inverse Gaussian pdf
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
