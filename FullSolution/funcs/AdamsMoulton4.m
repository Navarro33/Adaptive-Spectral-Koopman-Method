function x = AdamsMoulton4(f, x0, t, p)
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

dt = t(2) - t(1);
nt = numel(t);
nx = numel(x0);

x = zeros(nx, nt);
options = optimset('Display','off');

% first 5 steps
x(:, 1) = x0;

func = @(y) y - x(:, 1) - dt * f(t(2), y, p);
x(:, 2) = fsolve(func, x(:, 1), options);

func = @(y) y - x(:, 2) - dt * (1/2*f(t(3), y, p) + 1/2*f(t(2), x(:,2), p));
x(:, 3) = fsolve(func, x(:, 2), options);

func = @(y) y - x(:, 3) - dt * (5/12*f(t(4), y, p) + 8/12*f(t(3), x(:, 3), p) - 1/12*f(t(2), x(:, 2), p));
x(:, 4) = fsolve(func, x(:, 3), options);

func = @(y) y - x(:, 4) - dt * (9/24*f(t(5), y, p) + 19/24*f(t(4), x(:, 4), p) - 5/24*f(t(3), x(:, 3), p) + 1/24*f(t(2), x(:, 2), p));
x(:, 5) = fsolve(func, x(:, 4), options);

for s = 6:nt
    k2 = 646/720*f(t(s-1), x(:, s-1), p);
    k3 = -264/720*f(t(s-2), x(:, s-2), p);
    k4 = 106/720*f(t(s-3), x(:, s-3), p);
    k5 = -19/720*f(t(s-4), x(:, s-4), p);
    
    func = @(y) y - x(:, s-1) - dt * (251/720*f(t(s), y, p) + k2 + k3 + k4 + k5);
    x(:, s) = fsolve(func, x(:, s-1), options);
end