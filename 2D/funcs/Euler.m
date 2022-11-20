function x = Euler(f, x0, t, p)
% Euler method
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
x(:, 1) = x0;
for s = 1:(nt-1)
    x(:, s+1) = x(:, s) + dt * f(t(s), x(:,s), p);
end