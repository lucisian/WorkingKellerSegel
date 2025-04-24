function [Lap, Adv, ui, vi] = fnSetUpOperators(m, dx, dims)
    N = m^(dims);
    e = ones(m,1);
    Lap = spdiags([e, -2*e, e], [1, 0, -1], m, m);
    Adv = spdiags([e, -e], [1, -1], m, m);

    % Neumann BC
    Lap(1,1) = -1; Lap(end,end) = -1;
    Adv(1,1) = -1; Adv(end,end) = 1;

    ui = 1:N;
    vi = N+1:2*N;
end