function U0 = fnInitialCondition(N, a)
    rng(1);
    U0 = [ones(N,1); (1/a)*ones(N,1)] + 1e-2 * randn(2*N,1);
end
