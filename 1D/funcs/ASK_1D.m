function [xnT, n_decomp] = ASK_1D(dynamics, p, n, T, N, x0, r, fraction, option)
% The adaptive spectral Koopman method numerically solves autonomous
% dynamical systems, using the properties of the Koopman operator. This
% method also leverages ideas from the spectral collocation methods.
% This implemetation is based on one-dimensional systems. 
% 
% Args:
%       dynamics: the dynamics f of the system, a function handle
%       p: inherent parameters of the system
%       n: number of check points
%       T: terminal time
%       N: number of Gauss-Lobatto points 
%       x0: initial state
%       r: radius
%       fraction: fraction of radius, default 0.2
%       option: different type of polynomials for differentiation
%               1: Chebyshev (default) 2: Legendre
%
% Return:
%       xnT: numerical solution of the system at T
%       n_decomp: number of eigen-decompositions

%% parameter parsing
assert(isa(dynamics, 'function_handle'), 'TypeError: ');
assert(isstruct(p), 'TypeError: p must be a struct ...');
assert(n > 0, 'ValueError: n must be positive ...');
assert(T > 0, 'ValueError: T must be positive ...');
assert((N >= 3) && (mod(N, 2) == 1), 'ValueError: N must be at least 3 and odd ...');
assert(numel(x0) == numel(r), 'AttributeError: x0 and r must have the same length ...');
assert(numel(x0) == 1, 'AttributeError: x0 must be a scalar ...');

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

%% Algorithm
t = linspace(0, T, n+1);
lb = x0 - r;
ub = x0 + r;

% compute Gauss-Lobatto points and differentiation matrices
[x_GL_0, D_0] = compute_diffMat_1D(N, -1, 1, op);
[x_GL, D] = rescale_diffMat_1D(x_GL_0, D_0, lb, ub);

% compute finite dimensional Koopman approximation
[K, V, L] = approximate_Koopman_1D(N, x_GL, D, dynamics, p);

% compute projection coefficients
c = compute_coef_1D(x_GL, V);

% solve
xn = zeros(1, n+1);
xn(1) = x0;
n_decomp = 0;
t_ind = 1;
for j = 2:n+1
    % numerical solver
    x1 = numerical_solve_1D(t(j), t(t_ind), V, L, c);
    xn(j) = x1;
    
    % adaptive check
    flag = adaptive_check(x1, r, lb, ub, frac);
    
    % update
    if flag && (j < n+1)
        [lb, ub] = update_domain(x1, r);
        
        [x_GL, D] = rescale_diffMat_1D(x_GL_0, D_0, lb, ub);
        [K, V, L] = approximate_Koopman_1D(N, x_GL, D, dynamics, p);
        c = compute_coef_1D(x_GL, V);
        
        t_ind = j;
        n_decomp = n_decomp + 1;
    end
    
end

xnT = xn(end);
