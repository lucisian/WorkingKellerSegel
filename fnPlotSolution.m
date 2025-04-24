function fnPlotSolution(x, U, ui, dims)
    if dims == 1
        plot(x, U(end,ui), 'LineWidth', 3);
        xlabel('$x$','interpreter','latex');
        ylabel('$u(x,\textrm{end})$','interpreter','latex');
        title('Numerical Solution at End Time');
        set(gca, 'FontSize', 12);
        set(gca, 'TickLabelInterpreter', 'latex');
        axis tight
    else
        imagesc(reshape(U(end,ui), m, m));
        colorbar;
        axis off;
    end
end