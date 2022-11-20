function [F1, F2, F3] = KraichnanOrszagDynamics_slow(X, Y, Z, p, n_loop)

if nargin < 5
    n_loop = 5000;
end

F1 = [];
F2 = [];
F3 = [];
for k = 1:n_loop
    F1 = p.k1 * Y .* Z;
    F2 = p.k2 * X .* Z;
    F3 = p.k3 * (-2) * X .* Y;
end
