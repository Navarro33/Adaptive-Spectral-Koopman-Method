function [F1, F2] = SlowManifoldDynamics(X, Y, p)

F1 = p.alpha * X;
F2 = p.beta * (Y - X.^2);

