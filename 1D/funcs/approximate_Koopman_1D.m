function [K, V, L] = approximate_Koopman_1D(N, x_GL, D, dynamics, p)
% Finite dimensional approximation of the Koopman Operator and
% its eigen-decomposition for 1D systems. 
%
% Args:
%       N: number of Gauss-Lobatto points
%       x_GL: Gauss-Lobatto points
%       D: Differentiation matrix
%       dynamics: dynamics of the model, function handle
%       p: inherent parameters of the model, a struct
% Returns:
%       K: koopman approximation
%       V: eigenvector matrix
%       L: eigenvalue matrix

f = dynamics(x_GL, p);
K = diag(f) * D;

% Eigen-decomposition / spectral decomposition
[V, L] = eig(K);
