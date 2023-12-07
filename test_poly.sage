from sage.all import *

# 设置参数
q = 8380417
n = 256

# 创建多项式环
R.<x> = PolynomialRing(GF(q), implementation='NTL')

# 定义模
f = x**n + 1

# 创建商环
quotient_ring = QuotientRing(R, R.ideal(f))

# 创建两个多项式
poly1 = quotient_ring(sum(i * x**i for i in range(n)))
poly2 = quotient_ring(sum(((n - i) % q) * x**i for i in range(n)))

# 计算它们的和、差和积
sum_poly = poly1 + poly2
diff_poly = poly1 - poly2
prod_poly = poly1 * poly2


from random import randint
from param import *
from utils import *
from poly import polynomial

# 创建两个多项式
poly1 = polynomial([i % q for i in range(n)])
poly2 = polynomial([((n - i) % q) for i in range(n)])


# 计算它们的和、差和积
sum_poly1 = poly1 + poly2
diff_poly1 = poly1 - poly2
prod_poly1 = poly1 * poly2

print(sum_poly1.f == sum_poly.list())
print(diff_poly1.f == diff_poly.list())
print(prod_poly1.f == prod_poly.list())