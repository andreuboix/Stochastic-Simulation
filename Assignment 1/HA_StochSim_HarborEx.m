n = 90;
beta1 = 0.804590854705204;
beta2 =  0.832349335656409;
num_simulations = 1e6;
max_values = zeros(1, num_simulations);
mean_values = zeros(1, num_simulations);

for sim = 1:num_simulations
    W = zeros(1,90);
    T = weibull_rv(n, beta1);
    v = v_truncated(n);
    for i = 2:n
        if i == 2
            vu = v(1,i)/2;
            W(1,i) = max(0, vu - T(i));
        else
            if W(1,i-1) == 0
               vu = v(1,i)/2;
            else
               vu = v(1,i);
            end
            W(1,i) = max(0, W(1,i-1) + vu - T(i));
        end
    end
    
    max_values(sim) = max(W);
    mean_values(sim) = mean(W);
end
conf_int_max = quantile(max_values, [0.025 0.975]);
conf_int_mean = quantile(mean_values, [0.025 0.975]);

%fprintf('MAX is %f and AVERAGE is %f\n', max(W), mean(W))

fprintf('95%% Confidence Interval for MAX: [%f, %f]\n', conf_int_max(1), conf_int_max(2));
fprintf('95%% Confidence Interval for AVERAGE: [%f, %f]\n', conf_int_mean(1), conf_int_mean(2));

function X = weibull_rv(n, beta)
    U = rand(n, 1);          
    X = (-log(U)).^(1/beta);  
end


function v = v_truncated(n)
    v = randn(1,n);
    v(v > 1.5) = 1.5;
    v(v < 0.5) = 0.5;
end