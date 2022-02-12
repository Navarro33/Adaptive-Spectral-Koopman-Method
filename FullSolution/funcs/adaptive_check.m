function flag = adaptive_check(x_current, r, lb, ub, fraction)
% This function checks if the current solution falls out of the 
% acceptable domain. 
%
% Args:
%       x_current: current solution
%       r: radius
%       lb: lower bound
%       ub: upper bound
%       fraction: fraction of the radius
% Return:
%       flag: a bool 

flag = false;
if (x_current < lb + fraction*r) || (x_current > ub - fraction*r)
    flag = true;
end


