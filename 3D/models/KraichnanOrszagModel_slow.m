function g = KraichnanOrszagModel_slow(t, x, p, n_loop)
% Args:
%       t: time (here t is a dummy variable)
%       x: space var, 3D
%       p: parameters
% Returns:
%       g: derivative of x

if nargin < 4
    n_loop = 5000;
end

g = [];
for k = 1:n_loop
    g = [p.k1 * x(2) * x(3);
         p.k2 * x(1) * x(3);
         p.k3 * (-2) * x(1) * x(2)];
end
