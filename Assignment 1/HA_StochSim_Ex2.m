% Function to generate Weibull random variables using inversion method
function X = weibull_rv(n, beta)
    U = rand(n, 1);          % Uniform rv
    X = (-log(U)).^(1/beta);  % Apply inversion formula
end

% Theoretical Weibull PDF for comparison
function pdf_vals = weibull_pdf(x, beta)
    pdf_vals = beta * x.^(beta - 1) .* exp(-x.^beta);  % Weibull PDF
end

seed = 1234;
rng(seed);
% Parameters
beta = 1/2;      % Shape parameter (beta)
n = 1e7;       % Number of random variables to generate
bins = 3e2;       % Number of bins for the histogram

% Simulate Weibull random variables
X = weibull_rv(n, beta);

figure;
histogram(X,bins, 'Normalization','pdf');
hold off;
figure;
x_vals = linspace(0, max(X), 1000);
plot(x_vals, weibull_pdf(x_vals, beta), 'r-', 'LineWidth', 2);
hold off;

% Plot histogram of simulated values
figure;
histogram(X, bins, 'Normalization', 'pdf');
hold on;

% Line Plot of the theoretical density
x_vals = linspace(0, max(X), 1000);
plot(x_vals, weibull_pdf(x_vals, beta), 'r-', 'LineWidth', 2);

% Plot settings
title(['Histogram of Weibull rv (\beta = ', num2str(beta), ')']);
xlabel('x');
ylabel('Density');
legend('Simulated', 'Theoretical PDF');
grid on;
hold off;