function [xn, n_decomp] = ASK_Full_3D(dynamics, p, n, tspan, N, x0, r, fraction, option)
% The adaptive spectral Koopman method numerically solves autonomous
% dynamical systems, using the properties of the Koopman operator. This
% method also leverages ideas from the spectral collocation methods.
% This implemetation is based on three-dimensional systems, and solves for
% solutions at points of interest. 
% 
% Args:
%       dynamics: the dynamics f of the system, a function handle
%       p: inherent parameters of the system
%       n: number of check points
%       tspan: time span
%       N: number of Gauss-Lobatto points 
%       x0: initial state
%       r: radius
%       fraction: fraction of radius, default 0.2
%       option: different type of polynomials for differentiation
%               1: Chebyshev (default) 2: Legendre
%
% Return:
%       xn: numerical solutions of the system
%       n_decomp: number of eigen-decompositions

%% parameter parsing
assert(isa(dynamics, 'function_handle'), 'TypeError: ');
assert(isstruct(p), 'TypeError: p must be a struct ...');
assert(n > 0, 'ValueError: n must be positive ...');
assert((N >= 3) && (mod(N, 2) == 1), 'ValueError: N must be at least 3 and odd ...');
assert(numel(x0) == numel(r), 'AttributeError: x0 and r must have the same length ...');

x0 = reshape(x0, 3, 1);
frac = 0.2;
op = 1;

switch nargin
    case 8
        assert((fraction > 0) && (fraction <= 1), 'ValueError: fraction must be in (0, 1] ...');
        frac = fraction;
    case 9
        assert((fraction > 0) && (fraction <= 1), 'ValueError: fraction must be in (0, 1] ...');
        frac = fraction;
        assert((option == 1) || (option == 2), 'ValueError: option can be only 1 or 2 ...');
        op = option;
end

% early exit
if length(tspan) < 2
    xn = x0;
    n_decomp = 0;
    return;
end

%% algorithm
t = linspace(tspan(1), tspan(end), n+1);
lb = [x0(1)-r(1), x0(2)-r(2), x0(3)-r(3)];
ub = [x0(1)+r(1), x0(2)+r(2), x0(3)+r(3)];

% compute Gauss-Lobatto points and differentiation matrices
[x_GL_0, Ds_0] = compute_diffMat_3D(N, [-1, -1, -1], [1, 1, 1], op);
[x_GL, Ds] = rescale_diffMat_3D(x_GL_0, Ds_0, lb, ub);

% compute finite dimensional Koopman approximation
[K, V, L] = approximate_Koopman_3D(N, x_GL, Ds, dynamics, p);

% compute projection coefficients 
C = compute_coef_3D(x_GL, V);

% solver
xc = zeros(3, n+1);           % solutions at check points
xc = x0;
xn = zeros(3, length(tspan)); % solutions at time mesh points
xn(:,1) = x0;
n_decomp = 1;
t_ind = 1;                    % index of re-decomp check point 
tspan_ind = 1;                % index of time mesh point
for j = 2:n
    % numerical solver
    [x1, x2, x3] = numerical_solve_3D(t(j), t(t_ind), V, L, C);
    xc(:, j) = [x1, x2, x3];
    
    % adaptive check
    flag1 = adaptive_check(x1, r(1), lb(1), ub(1), frac);
    flag2 = adaptive_check(x2, r(2), lb(2), ub(2), frac);
    flag3 = adaptive_check(x3, r(3), lb(3), ub(3), frac);
    
    % update
     if (flag1 || flag2 || flag3)
%         % solve at time mesh points
%         tspan_ind_new = find((tspan < t(j)+1e-15), 1, 'last');
%         if tspan_ind_new > tspan_ind
%             for i = tspan_ind+1 : tspan_ind_new
%                 [sol1, sol2, sol3] = numerical_solve_3D(tspan(i), t(t_ind), V, L, C);
%                 xn(:, i) = [sol1, sol2, sol3];
%             end
%             tspan_ind = tspan_ind_new;
%         end
        
        tspan_ind_new = find((tspan(tspan_ind+1:end) < t(j)+1e-15), 1, 'last') + tspan_ind;
        if ~isempty(tspan_ind_new)
            ct = tspan(tspan_ind+1:tspan_ind_new);
            [sol1, sol2, sol3] = numerical_solve_3D(ct, t(t_ind), V, L, C);
            xn(:, tspan_ind+1:tspan_ind_new) = [sol1; sol2; sol3];
            
            tspan_ind = tspan_ind_new;
        end
        
        [lb(1), ub(1)] = update_domain(x1, r(1));
        [lb(2), ub(2)] = update_domain(x2, r(2));
        [lb(3), ub(3)] = update_domain(x3, r(3));
        
        [x_GL, Ds] = rescale_diffMat_3D(x_GL_0, Ds_0, lb, ub);
        [K, V, L] = approximate_Koopman_3D(N, x_GL, Ds, dynamics, p);
        C = compute_coef_3D(x_GL, V);
        
        t_ind = j;
        n_decomp = n_decomp + 1;
    end
        
end

% solve at the rest of the time mesh points
% if tspan_ind < length(tspan)
%     for i = tspan_ind+1 : length(tspan)
%         [sol1, sol2, sol3] = numerical_solve_3D(tspan(i), t(t_ind), V, L, C);
%         xn(:, i) = [sol1, sol2, sol3]';
%     end
% end

if tspan_ind < length(tspan)
    ct = tspan(tspan_ind+1:length(tspan));
    [sol1, sol2, sol3] = numerical_solve_3D(ct, t(t_ind), V, L, C);
    xn(:, tspan_ind+1:length(tspan)) = [sol1; sol2; sol3];
end

