function g = KraichnanOrszagModel(t, x, p)
% Args:
%       t: time (here t is a dummy variable)
%       x: space var, 3D
%       p: parameters
% Returns:
%       g: derivative of x

g = [p.k1 * x(2) * x(3);
     p.k2 * x(1) * x(3);
     p.k3 * (-2) * x(1) * x(2)];
