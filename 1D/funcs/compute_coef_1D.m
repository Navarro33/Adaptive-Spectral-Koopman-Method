function c = compute_coef_1D(x_GL, V)
% Compute the coefficients of the linear combination of eigenvectors. 
% This linear combination projects the obvervable into the vector space 
% spanned by eigenvectors.
% 
% Args:
%       x_GL: Gauss-Lobatto points
%       V: matrix of eigenvectors
% Returns:
%       c: coefficient vector

c = V \ x_GL;
