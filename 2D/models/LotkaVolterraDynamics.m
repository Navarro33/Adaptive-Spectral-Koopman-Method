function [F1, F2] = LotkaVolterraDynamics(X, Y, p)

F1 = p.alpha * X - p.beta * X .* Y;
F2 = p.delta * X .* Y - p.gamma * Y;

