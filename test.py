import sympy as sp

ex = sp.Symbol('ex')
ey = sp.Symbol('ey')
ez = sp.Symbol('ez')

x = sp.Symbol('x')
y = sp.Symbol('y')
z = sp.Symbol('z')

expr = (x*ex+y*ey+z*ez)*(x*ex+y*ey+z*ez)
expr = sp.expand(expr).subs({ex*ey:sp.cos(x-y), ex*ez:sp.cos(z-x), ey*ez:sp.cos(y-z), ex**2:1, ey**2:1, ez**2:1})
print(expr)