function C = compute_coef_2D(x_GL, V)
% Compute the coefficients of the linear combination of eigenmatrices. 
% This linear combination projects the obvervable into the vector space 
% spanned by eigenmatrices.
% 
% Args:
%       x_GL: Gauss-Lobatto points
%       V: matrix of eigenvectors
% Returns:
%       C: coefficient matrix whose columns are the coefficients

[m, ~] = size(x_GL);

x1_GL = x_GL(:, 1);
x2_GL = x_GL(:, 2);

% c1 = V \ repmat(x1_GL, m, 1);
% c2 = V \ repelem(x2_GL, m);
tmp1 = repmat(x1_GL, m, 1);
tmp2 = repelem(x2_GL, m);

C = V \ [tmp1, tmp2];
