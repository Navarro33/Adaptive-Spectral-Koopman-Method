function g = LimitCycle(t, x, p)
% Compute the derivatives of the limit cycle.
% Args:
%       t: time (here t is a dummy variable)
%       x: space var, 2D
%       p: parameters
% Returns:
%       g: derivative of x

g = [- x(1) - x(2) + p.k * x(1) / sqrt(x(1)^2 + x(2)^2);
       x(1) - x(2) + p.k * x(2) / sqrt(x(1)^2 + x(2)^2)];


