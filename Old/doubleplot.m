[X,Y] = meshgrid(-2:.1:2,-2:.1:2); % Generate domain.
Z = 0.053772*X+0.570473*X.^(3)-0.394611*X.^(5)+(1.07049*X-1.48409*X.^(3))*Y.^(2)-290.144*X.^(4)*Y; % Find function value everywhere in the domain.

contour(X,Y,Z,[0 0]) % Plot the isoline where the function value is 4.