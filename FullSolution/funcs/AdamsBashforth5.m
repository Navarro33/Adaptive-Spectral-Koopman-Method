function x = AdamsBashforth5(f, x0, t, p)
% 5-step Adams-Bashforth method
%
% Args:
% Args:
%       f: derivatives, a function of time and state
%       x0: initial condition
%       t: time mesh
%       p: inherent parameters of f
% Returns:
%       x: numerical solution

dt = t(2) - t(1);
nt = numel(t);
nx = numel(x0);

x = zeros(nx, nt);

% first 5 steps
x(:, 1) = x0;
x(:, 2) = x(:, 1) + dt * f(t(1), x(:, 1), p);
x(:, 3) = x(:, 2) + dt * (3/2*f(t(2), x(:, 2), p) - 1/2*f(t(1), x(:, 1), p));
x(:, 4) = x(:, 3) + dt * (23/12*f(t(3), x(:, 3), p) - 16/12*f(t(2), x(:, 2), p) + 5/12*f(t(1), x(:, 1), p));
x(:, 5) = x(:, 4) + dt * (55/24*f(t(4), x(:, 4), p) - 59/24*f(t(3), x(:, 3), p) + 37/24*f(t(2), x(:, 2), p) - 9/24*f(t(1), x(:, 1), p));

for s = 6:nt
    k1 = 1901/720*f(t(s-1), x(:, s-1), p);
    k2 = -2774/720*f(t(s-2), x(:, s-2), p);
    k3 = 2616/720*f(t(s-3), x(:, s-3), p);
    k4 = -1274/720*f(t(s-4), x(:, s-4), p);
    k5 = 251/720*f(t(s-5), x(:, s-5), p);
    x(:, s) = x(:, s-1) + dt * (k1+k2+k3+k4+k5);
end