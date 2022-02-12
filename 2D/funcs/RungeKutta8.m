function x = RungeKutta8(f, x0, t, p)
% Fehlberg's 8th Order Runge-Kutta methods.
%
% Args:
%       f: derivatives, a function of time and state
%       x0: initial condition
%       t: time mesh
%       p: inherent parameters of f
% Returns:
%       x: numerical solution

h = t(2) - t(1);
nt = numel(t);
nx = numel(x0);

x = zeros(nx, nt);
x(:, 1) = x0;
for s = 1:(nt - 1)
    ts = t(s);
    xs = x(:,s);
    
    k1 = f(ts, xs, p);
    k2 = f(ts + 2*h/27, xs + 2*h*k1/27, p);
    k3 = f(ts + h/9, xs + h/36*(k1 + 3*k2), p);
    k4 = f(ts + h/6, xs + h/24*(k1 + 3*k3), p);
    k5 = f(ts + 5*h/12, xs + h/48*(20*k1 - 75*k3 + 75*k4), p);
    k6 = f(ts + h/2, xs + h/20*(k1 + 5*k4 + 4*k5), p);
    k7 = f(ts + 5*h/6, xs + h/108*(-25*k1 + 125*k4 - 260*k5 + 250*k6), p);
    k8 = f(ts + h/6, xs + h*(31/300*k1 + 61/225*k5 - 2/9*k6 + 13/900*k7), p);
    k9 = f(ts + 2*h/3, xs + h*(2*k1 - 53/6*k4 + 704/45*k5 - 107/9*k6 + 67/90*k7 + 3*k8), p);
    k10 = f(ts + h/3, xs + h*(-91/108*k1 + 23/108*k4 - 976/135*k5 + 311/54*k6 - 19/60*k7 + 17/6*k8 - 1/12*k9), p);
    k11 = f(ts + h, xs + h*(2383/4100*k1 - 341/164*k4 + 4496/1025*k5 - 301/82*k6 + 2133/4100*k7 + 45/82*k8 + 45/164*k9 + 18/41*k10), p);
%     k12 = f(ts, xs + h*(3/205*k1 - 6/41*k6 - 3/205*k7 - 3/41*k8 + 3/41*k9 + 6/41*k10), p);
%     k13 = f(ts + h, xs + h*(-1777/4100*k1 - 341/164*k4 + 4496/1025*k5 - 289/82*k6 + 2193/4100*k7 + 51/82*k8 + 33/164*k9 + 12/41*k10 + k12), p);
    
    dx = h * (41/840*k1 + 34/105*k6 + 9/35*k7 + 9/35*k8 + 9/280*k9 + 9/280*k10 + 41/840*k11);
    x(:, s+1) = xs + dx;
end