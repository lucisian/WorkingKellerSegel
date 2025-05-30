function fnSpaceTimePlot(U, x, T)
% SpaceTimePlot - Plots a space-time diagram for u(x,t)
%
% Inputs:
%   U - solution matrix from ode15s, size (length(T) x 2*N)
%   x - spatial grid (vector of length m)
%   T - time vector (same length as time steps in U)
%
% Assumes u is the first half of the state vector.

    [~, ~] = meshgrid(x, T);  % (not strictly necessary for imagesc)

    figure;
    imagesc(T, x, U(:, 1:end/2)');  % Plot u(x,t)
    set(gca, 'YDir', 'normal');    % So x increases from bottom to top
    colorbar;

    xlabel('$t$','interpreter','latex');
    ylabel('$x$','interpreter','latex');
    title('Space-Time Plot');
    colormap(parula);              % Default MATLAB colormap
    set(gca, 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex');
end