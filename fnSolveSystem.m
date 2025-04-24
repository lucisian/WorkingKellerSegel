function [T,U] = fnSolveSystem(F, U0, Lap, N, T)
    JacSparse = sparse([Lap, Lap; eye(N), Lap]);
    odeOptions = odeset('JPattern', JacSparse);
    [T,U] = ode15s(F, T, U0, odeOptions);
end