function xn = numerical_solve_1D(ct, lt, V, L, c)
% Numerically solve the equations using the eigenvectors 
% and eigenvalues given by the Koopman approx. 
%
% In particular, 
%   x = sum_i c_i * v_i(x_init) * exp(lambda* t),
% where v_i is eigenfunction value given by the eigenvector.
%
% Args:
%       ct: current time
%       lt: last time (of eigen-decomposition)
%       V: matrix of eigenvectors
%       L: matrix of eigenvalues
%       c: coefficients
% Returns:
%       xn: numerical solution of the first component



[nrow, ~] = size(V);
idx = (nrow + 1) / 2;      % middle point

% solve a linear dynamics
lambda = diag(L);
v = V(idx, :);
delta_t = ct - lt;

xn = real( v * (c .* exp(lambda * delta_t)) );





