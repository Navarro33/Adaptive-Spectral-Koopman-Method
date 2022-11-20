function [F1, F2] = SimplePendulumDynamics_slow(X, Y, p)

F1 = Y;
F2 = - (p.g / p.L) * sin_taylor(X);
