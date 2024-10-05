function Convergence
    seed = 1234;
    rng(seed);
    n = 1e6;
    tol = 1e-6;
    p = pi;

    [n, result] = HA_StochSim_Ex1(n);
    while (abs(result - p) > tol)
        n = n + 1;
        [n, result] = HA_StochSim_Ex1(n);
        %disp([n, result]);
    end
    X = sprintf("The number of samples is %0.5e and the result to which it approaches is %0.7d", n, result);
    disp(X);
end
