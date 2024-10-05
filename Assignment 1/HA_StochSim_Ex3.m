function HA_StochSim_Ex3
    % Parameters of the Inverse Gaussian distribution
    c = 1; 
    xi = 1;
    n_samples = 1e6; % Number of samples wanted  

    lambda = xi^2/2;  % Parameter for the exponential distribution
    M = 2.4944; % Tuned scaling constant
    
    % Initialize sample array
    samples = zeros(n_samples, 1);
    num_accepted = 0;
    
    % Generate samples using the acceptance-rejection method
    while num_accepted < n_samples
        X_star = exprnd(1 / lambda);
        
        U = rand;
        
        accept_prob = inverse_gaussian_pdf(X_star, c, xi) / (M * exponential_proposal(X_star, lambda));
        
        % Accept or reject the candidate
        if U <= accept_prob
            num_accepted = num_accepted + 1;
            samples(num_accepted) = X_star;
        end
    end
    
    % Compute the sample mean, squared mean, and variance
    sample_mean = mean(samples);
    sample_squared_mean = mean(samples.^2);
    sample_variance = var(samples);
    
    fprintf('Sample Mean: %.4f\n', sample_mean);
    fprintf('Sample Squared Mean: %.4f\n', sample_squared_mean);
    fprintf('Sample Variance: %.4f\n', sample_variance);
    
    % Compare with theoretical values
    theoretical_mean = c / xi;
    theoretical_variance = c / xi^3;
    
    fprintf('Theoretical Mean: %.4f\n', theoretical_mean);
    fprintf('Theoretical Squared Mean: %.4f\n', theoretical_mean^2 + theoretical_variance);
    fprintf('Theoretical Variance: %.4f\n', theoretical_variance);
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
