function C = compute_coef_3D(x_GL, V)
% Compute the coefficients of the linear combination of eigentensors. 
% This linear combination projects the obvervable into the vector space 
% spanned by eigentensors.
% 
% Args:
%       x_GL: Gauss-Lobatto points
%       V: matrix of eigenvectors
% Returns:
%       C: coefficient matrix whose columns are the coefficients

[m, ~] = size(x_GL);

x1_GL = x_GL(:, 1);
x2_GL = x_GL(:, 2);
x3_GL = x_GL(:, 3);

% c1 = V \ repmat(x1_GL, m^2, 1);
% c2 = V \ repelem(repmat(x2_GL, m, 1), m);
% c3 = V \ repelem(x3_GL, m^2);
tmp1 = repmat(x1_GL, m^2, 1);
tmp2 = repelem(repmat(x2_GL, m, 1), m);
tmp3 = repelem(x3_GL, m^2);

C = V \ [tmp1, tmp2, tmp3];
