function x = AdamsMoulton4(f, x0, t, p, warm_up)
% 4-step Adams-Moulton method
%
% Args:
% Args:
%       f: derivatives, a function of time and state
%       x0: initial condition
%       t: time mesh
%       p: inherent parameters of f
% Returns:
%       x: numerical solution

if nargin < 5
    warm_up = true;
end

dt = t(2) - t(1);
nt = numel(t);
nx = numel(x0);

x = zeros(nx, nt);
options = optimset('Display','off', 'TolFun', 1e-15);

% first 5 steps
if warm_up
    % use Runge-Kutta 4 to solve first steps
    x(:, 1) = x0;
    for s = 1:4
        k1 = dt * f(t(s), x(:, s), p);
        k2 = dt * f(t(s) + dt/2, x(:, s) + k1/2, p);
        k3 = dt * f(t(s) + dt/2, x(:, s) + k2/2, p);
        k4 = dt * f(t(s) + dt, x(:, s) + k3, p);

        dx = (k1 + 2*k2 + 2*k3 + k4) / 6;
        x(:, s+1) = x(:, s) + dx;
    end
else
    x(:, 1) = x0;

    func = @(y) y - x(:, 1) - dt * f(t(2), y, p);
    x(:, 2) = fsolve(func, x(:, 1), options);

    func = @(y) y - x(:, 2) - dt/2 * (1*f(t(3), y, p) + 1*f(t(2), x(:,2), p));
    x(:, 3) = fsolve(func, x(:, 2), options);

    func = @(y) y - x(:, 3) - dt/12 * (5*f(t(4), y, p) + 8*f(t(3), x(:, 3), p) - 1*f(t(2), x(:, 2), p));
    x(:, 4) = fsolve(func, x(:, 3), options);

    func = @(y) y - x(:, 4) - dt/24 * (9*f(t(5), y, p) + 19*f(t(4), x(:, 4), p) - 5*f(t(3), x(:, 3), p) + 1*f(t(2), x(:, 2), p));
    x(:, 5) = fsolve(func, x(:, 4), options);
end

for s = 6:nt
    k2 = 646*f(t(s-1), x(:, s-1), p);
    k3 = -264*f(t(s-2), x(:, s-2), p);
    k4 = 106*f(t(s-3), x(:, s-3), p);
    k5 = -19*f(t(s-4), x(:, s-4), p);
    
    func = @(y) y - x(:, s-1) - dt/720 * (251*f(t(s), y, p) + k2 + k3 + k4 + k5);
    x(:, s) = fsolve(func, x(:, s-1), options);
end