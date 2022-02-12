function [x_GL_new, D_new] = rescale_diffMat_1D(x_GL, D, lb, ub)
% This function rescales the Gauss-Lobatto points and the differentiation
% matrices from [-1, 1] to [lb, ub].
%
% Args:
%       x_GL: Gauss-Lobatto points
%       D: Differentiation matix
%       lb: lower bound
%       ub: upper bound
% Return:
%       x_GL_new: Rescaled Gauss-Lobatto points
%       D_new: Rescaled differentiation matix


x_GL_new = (ub - lb)/2 * (x_GL + 1) + lb;
D_new = 2/(ub - lb) * D;
