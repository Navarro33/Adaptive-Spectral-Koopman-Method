function g = SlowManifold(t, x, p)
% Compute the derivatives of the slow manifold.
% Args:
%       t: time (here t is a dummy variable)
%       x: space var, 2D
%       p: parameters
% Returns:
%       g: derivative of x

g = [p.alpha * x(1);
     p.beta * ( x(2) - x(1)^2 )];
