from sage.all import *

q = 8380417
n = 256

R.<x> = PolynomialRing(GF(q), implementation='NTL')
f = x**n + 1
quotient_ring = QuotientRing(R, R.ideal(f))

poly1 = quotient_ring(sum(i * x**i for i in range(n)))
poly2 = quotient_ring(sum(((n - i) % q) * x**i for i in range(n)))

sum_poly = poly1 + poly2
diff_poly = poly1 - poly2
prod_poly = poly1 * poly2


from random import randint
from param import *
from utils import *
from poly import polynomial

poly1 = polynomial([i % q for i in range(n)])
poly2 = polynomial([((n - i) % q) for i in range(n)])

sum_poly1 = poly1 + poly2
diff_poly1 = poly1 - poly2
prod_poly1 = poly1 * poly2

print(sum_poly1.f == sum_poly.list())
print(diff_poly1.f == diff_poly.list())
print(prod_poly1.f == prod_poly.list())
