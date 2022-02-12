function x = RungeKutta4(f, x0, t, p)
% 4th Order Runge-Kutta methods.
%
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
x(:, 1) = x0;
for s = 1:(nt-1)
    k1 = dt * f(t(s), x(:, s), p);
    k2 = dt * f(t(s) + dt/2, x(:, s) + k1/2, p);
    k3 = dt * f(t(s) + dt/2, x(:, s) + k2/2, p);
    k4 = dt * f(t(s) + dt, x(:, s) + k3, p);
    
    dx = (k1 + 2*k2 + 2*k3 + k4) / 6;
    x(:, s+1) = x(:, s) + dx;
end

