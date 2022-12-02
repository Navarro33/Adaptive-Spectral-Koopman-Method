function y = solve_limit_cycle(t, x, p) 
assert(p.k == 1, 'This function only applies to p.k = 1 case');
t = reshape(t, [1, length(t)]);
x = reshape(x, [], 2);

y_1 = (1 - (1 - sqrt(x(:,1).^2 + x(:,2).^2))*exp(-t) ).*cos(t + atan(x(:,2)./x(:,1)));
y_2 = (1 - (1 - sqrt(x(:,1).^2 + x(:,2).^2))*exp(-t) ).*sin(t + atan(x(:,2)./x(:,1)));

if size(t,2) == 1
    y = [y_1, y_2];
elseif size(t, 2) > 1 && size(x, 1) == 1
    y = [y_1; y_2];
else
    y = cat(3, y_1, y_2);
end


