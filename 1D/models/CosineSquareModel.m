function g = CosineSquareModel(t, x, p)
% Compute the derivatives of the cosine model.
% Args:
%       t: time (here t is a dummy variable)
%       x: space var, 2D
%       p: parameters, a struct with alpha, beta, gamma, delta
% Returns:
%       g: derivative of x
g = p.k * cos(x)^2;