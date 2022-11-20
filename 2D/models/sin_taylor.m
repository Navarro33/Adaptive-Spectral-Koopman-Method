function y = sin_taylor(x, N)
% A slow implementation of Taylor expansion of sin(x) to 2N terms

if nargin < 2
    N = 200;
end

y = zeros(size(x));
for k = 0:N
    tmp = (x.^(2*k+1) ./ factorial(2*k+1));
    if isnan(tmp)
        tmp = 0;
    end
    y = y + (-1)^k * tmp;
end
