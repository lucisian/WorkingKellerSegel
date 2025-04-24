clear;
% Define fixed and varying parameters
a = 1; D = 1; p2 = 1; p4 = 0; epsilon = 0.2;
run('parameters5.m');

% Setup operators
[Lap, Adv, ui, vi] = fnSetUpOperators(m, dx, dims);

% Compute p, define reaction function
F = fnDefineKinetics(a, D, p, Lap, Adv, dx, dims, ui, vi);

% Initial condition
U0 = fnInitialCondition(N, a);

% Solve PDE
[T, U] = fnSolveSystem(F, U0, Lap, N, T);

% Plot results
fnPlotSolution(x, U, ui, dims);
fnSpaceTimePlot(U,x,T);
%fnPlotKymograph(U, x, T, ui);
fnPlotDominantFrequency(U, ui, T, m);