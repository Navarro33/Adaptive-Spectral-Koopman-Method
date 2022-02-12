function g = LotkaVolterraModel(t, x, p)
% Compute the derivatives of the Lotka-Volterra model.
% Args:
%       t: time (here t is a dummy variable)
%       x: space var, 2D
%       p: parameters
% Returns:
%       g: derivative of x
g = [p.alpha*x(1) - p.beta*x(1)*x(2);
     p.delta*x(1)*x(2) - p.gamma*x(2)];
    
