import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy import Derivative

t = sp.Symbol('t')
x0 = sp.Function('x0')(t)
x1 = x0.diff(t)
x2 = x1.diff(t)

X1 = sp.Function('X1')(t)
X2 = sp.Symbol('X2')


y0 = sp.Function('y0')(t)
y1 = y0.diff(t)
y2 = y1.diff(t)

Y1 = sp.Function('Y1')(t)
Y2 = sp.Symbol('Y2')

g = sp.Symbol('g')
m_1 = sp.Symbol('m_1')
m_2 = sp.Symbol('m_2')
l_1 = sp.Symbol('l_1')
l_2 = sp.Symbol('l_2')

E_k1 = m_1 * l_1 * l_1 * x1 * x1
E_k2 = m_2 * (x1 * x1 * l_1 * l_1 + y1 * y1 * l_2 * l_2 + 2 * x1 * y1 * l_2 * l_2 * sp.cos(y0 - x0))

E_p1 = m_1 * g * l_1 * (1-sp.cos(x0))
E_p2 = m_2 * g * (l_1 * (1 - sp.cos(x0)) + l_2 * (1 - sp.cos(y0)))

T = E_k1 + E_k2
V = E_p1 + E_p2
L = T - V

a = (L.diff(x1).diff(t) - L.diff(x0))
print(a)
a = a.subs(x0.diff(t),X1).subs(y0.diff(t),Y1)
a = a.subs(X1.diff(t),X2).subs(Y1.diff(t),Y2)
print(a)
b = (L.diff(y1).diff(t) - L.diff(y0))
print(b)
b = b.subs(x0.diff(t),X1).subs(y0.diff(t),Y1)
b = b.subs(X1.diff(t),X2).subs(Y1.diff(t),Y2)
print(sp.simplify(b))


print(sp.solve([sp.Eq(a,0), sp.Eq(b,0)],X2,Y2))
