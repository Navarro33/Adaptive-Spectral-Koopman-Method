function f = CosineSquareDynamics_slow(x, p)

y = cos2_taylor(x);
f = p.k * y;
