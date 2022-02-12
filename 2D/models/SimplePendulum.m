function g = SimplePendulum(t, x, p)
% Compute the derivatives of the simple pendulum.
% xtt = - g/L sin(x) --> x1t = x2 and x2t = - g/L sin(x1)
%
% Define x = (x1, x2) are two-dimensional and t is a dummy.
% Returns:
%           g: derivative of x

g = [ x(2);
      - (p.g / p.L) * sin(x(1))];

