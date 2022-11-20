function x = AdamsBashforth5(f, x0, t, p, warm_up)
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

if nargin < 5
    warm_up = true;
end

dt = t(2) - t(1);
nt = numel(t);
nx = numel(x0);

x = zeros(nx, nt);

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
%     for s = 1:4
%         ts = t(s);
%         xs = x(:,s);
% 
%         k1 = f(ts, xs, p);
%         k2 = f(ts + dt/4, xs + dt*k1/4, p);
%         k3 = f(ts + 3*dt/8, xs + dt*(3*k1/32 + 9*k2/32), p);
%         k4 = f(ts + 12*dt/13, xs + dt*(1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197), p);
%         k5 = f(ts + dt, xs + dt*(439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104), p);
%         %k6 = f(ts + dt/2, xs + dt*(-8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40), p);
% 
%         dx = dt * (25/216 * k1 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5);
%         x(:, s+1) = xs + dx;
%     end
else
    % original first steps of 5-step Adams-Bashforth
    x(:, 1) = x0;
    x(:, 2) = x(:, 1) + dt * f(t(1), x(:, 1), p);
    x(:, 3) = x(:, 2) + dt/2 * (3*f(t(2), x(:, 2), p) - 1*f(t(1), x(:, 1), p));
    x(:, 4) = x(:, 3) + dt/12 * (23*f(t(3), x(:, 3), p) - 16*f(t(2), x(:, 2), p) + 5*f(t(1), x(:, 1), p));
    x(:, 5) = x(:, 4) + dt/24 * (55*f(t(4), x(:, 4), p) - 59*f(t(3), x(:, 3), p) + 37*f(t(2), x(:, 2), p) - 9*f(t(1), x(:, 1), p));
end

for s = 6:nt
    k1 = 1901*f(t(s-1), x(:, s-1), p);
    k2 = -2774*f(t(s-2), x(:, s-2), p);
    k3 = 2616*f(t(s-3), x(:, s-3), p);
    k4 = -1274*f(t(s-4), x(:, s-4), p);
    k5 = 251*f(t(s-5), x(:, s-5), p);
    x(:, s) = x(:, s-1) + dt/720 * (k1+k2+k3+k4+k5);
end