function F = fnDefineKinetics(a, D, p, Lap, Adv, dx, dims, ui, vi)
    f = @(u,v) u.*(1-u);
    g = @(u,v) u - a*v;

    if dims == 1
        Lap = (1/dx^2)*Lap;
        Adv = (1/(2*dx))*Adv;
        chi = @(U) U ./ (1 + U.^2);
        F = @(t,U)[...
            f(U(ui),U(vi)) + Lap*U(ui) - ...
            p*(chi(U(ui)).*(Lap*U(vi)) + (Adv*(chi(U(ui)))).*(Adv*U(vi))); ...
            g(U(ui),U(vi)) + D*Lap*U(vi)];
    else
        error("2D not yet implemented here.");
    end
end