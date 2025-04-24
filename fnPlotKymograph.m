function fnPlotKymograph(U, x, T, ui)
% Set up for frequency tracking
%m=4000;
fs = m; % sampling frequency in space
Nl = m; % number of spatial points
freq = 0:fs/Nl:fs/2; % frequency axis (positive freqs only)
dominant_freq = zeros(length(T),1); % to store dominant frequency at each time

% Loop through each time step
for i = 1:length(T)
    u = U(i, ui); % extract u(x) at current time
    u = u - mean(u); % zero-mean for FFT
    
    xdft = fft(u);
    xdft = xdft(1:Nl/2+1); % one-sided FFT
    psdx = (1/(fs*Nl)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    
    % Find dominant frequency (highest power)
    [~, idx] = max(psdx);
    dominant_freq(i) = freq(idx);
end

% Plot the dominant frequency over time
figure;
plot(T, 2*dominant_freq, 'LineWidth', 2);
xlabel('Time');
ylabel('Dominant Frequency');
title('Dominant Spatial Frequency Over Time');
grid on;
% PlotKymograph - Plots a space-time diagram (kymograph) for variable u(x,t)
%
% Usage:
%   PlotKymograph(U, x, T, ui)
%
% Inputs:
%   U  - Solution matrix from ode15s, size (length(T) × 2*N)
%   x  - Spatial domain vector
%   T  - Time points vector
%   ui - Index vector corresponding to u in the solution

    % Create the space-time image plot
   % imagesc(T, x, U(:, ui)'); % transpose so space is vertical, time is horizontal
   % set(gca, 'YDir', 'normal'); % ensure y-axis is not flipped

    % Axis labels with LaTeX interpreter
   % xlabel('$t$', 'Interpreter', 'latex');
   % ylabel('$x$', 'Interpreter', 'latex');

    % Set color map and shading
   % colormap(viridis);  % requires viridis colormap in path, or replace with 'parula'
   % shading interp;

    % Colorbar settings
   % c = colorbar;
   % c.Label.String = '$u$';
   % c.TickLabelInterpreter = 'latex';
   % c.Label.Interpreter = 'latex';

    % Style settings
   % set(gca, 'TickLabelInterpreter', 'latex');
   % set(gca, 'FontSize', 24);

   % title('Kymograph of $u(x,t)$', 'Interpreter', 'latex');
end