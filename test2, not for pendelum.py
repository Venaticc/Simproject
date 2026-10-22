import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy import Derivative

t = sp.Symbol('t')
x0 = sp.Function('x0')(t)
x1 = x0.diff(t)
x2 = x1.diff(t)

X1 = sp.Symbol('X1')
X2 = sp.Symbol('X2')


y0 = sp.Function('y0')(t)
y1 = y0.diff(t)
y2 = y1.diff(t)

Y1 = sp.Symbol('Y1')
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

A = L.diff(x1).diff(t) - L.diff(x0)
print(sp.simplify(A))
B = L.diff(y1).diff(t) - L.diff(y0)
print(sp.simplify(B))

print(sp.solve([sp.Eq(A,0), sp.Eq(B,0)],x2,y2))