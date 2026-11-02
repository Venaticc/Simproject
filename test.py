import sympy as sp

x = sp.Symbol('x')
y = sp.Symbol('y')
z = sp.Symbol('z')

ex = sp.Symbol('ex')
ey = sp.Symbol('ey')
ez = sp.Symbol('ez')

vectordict = {ex*ey:sp.cos(x-y), ex*ez:sp.cos(z-x), ey*ez:sp.cos(y-z), ex**2:1, ey**2:1, ez**2:1}

l_x = sp.Symbol('l_x')
l_y = sp.Symbol('l_y')
l_z = sp.Symbol('l_z')
m_x = sp.Symbol('m_x')
m_y = sp.Symbol('m_y')
m_z = sp.Symbol('m_z')

symvaldict = {l_x:1, l_y:1, l_z:1, m_x:1, m_y:1, m_z:1}

expr = (x*ex+y*ey+z*ez)*(x*ex+y*ey+z*ez)
expr = sp.expand(expr).subs(vectordict)
print(expr)