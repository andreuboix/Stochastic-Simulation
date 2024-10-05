function [n,result] = HA_StochSim(n)
    % We define U1 as the n-vector for which, at every position,
    % there is a sample of a uniform distribution on [-1, 1].
    U1 = -1 + 2 * rand(n,1);
    
    % We define U2 as the n-vector for which, at every position,
    % there is a sample of a uniform distribution on [0, 1].
    U2 = rand(n,1);
    
    U_samples = [U1, U2];
    
    sum_value = 0;
    
    % Loop over each sample and check if it lies within the unit circle
    for i = 1:n
        U1_i = U_samples(i, 1);
        U2_i = U_samples(i, 2);
        
        % Check the condition (U1_i)^2 + (U2_i)^2 <= 1
        if (U1_i^2 + U2_i^2) <= 1
            sum_value = sum_value + 1;
        end
    end
    
    % Compute the final result
    result = (4 / n) * sum_value;
    
    %disp('Number of samples and computed result:');
    %disp([n, result]);
end