function [lb_new, ub_new] = update_domain(x_current, r)
% Updates the domain for computing Gauss-Lobatto points and the
% differentiation matrix.
%
% Args:
%       x_current: current solution
%       r: radius
% Return:
%       lb_new: updated lower bound
%       ub_new: updated upper bound

lb_new = x_current - r;
ub_new = x_current + r;


