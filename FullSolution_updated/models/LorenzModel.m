function g = LorenzModel(t, x, p)
% Args:
%       t: time (here t is a dummy variable)
%       x: space var, 3D
%       p: parameters
% Returns:
%       g: derivative of x

g = [p.alpha * (x(2) - x(1));
     x(1) * (p.beta - x(3)) - x(2);
     x(1) * x(2) - p.gamma * x(3)];
