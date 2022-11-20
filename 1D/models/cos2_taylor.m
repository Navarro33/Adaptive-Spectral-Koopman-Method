function y = cos2_taylor(x, N)
% 
% A slow implementation of Taylor expansion of cos(x)^2 to 2N terms

if nargin < 2
    N = 500;
end
y = ones(size(x));
for k = 1:N
    tmp = x.^(2*k) * (2^(-1+2*k) / factorial(2*k));
    if isnan(tmp)
        tmp = 0;
    end
    y = y + (-1)^k * tmp;
end

