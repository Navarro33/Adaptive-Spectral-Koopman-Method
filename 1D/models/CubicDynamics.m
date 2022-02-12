function f = CubicDynamics(x, p)

f = p.k ./ (3 .* x.^2);

