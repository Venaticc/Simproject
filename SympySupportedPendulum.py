import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy import Derivative

t = sp.Symbol('t') # time

x0 = sp.Function('x0')(t)   # theta
x1 = x0.diff(t)             # theta prime
x2 = x1.diff(t)             # theta bis

X1 = sp.Function('X1')(t)   # replacement function, magic trick
X2 = sp.Symbol('X2')

y0 = sp.Function('y0')(t)   # phi
y1 = y0.diff(t)             # phi prime
y2 = y1.diff(t)             # phi bis

Y1 = sp.Function('Y1')(t)   # also magic tricks
Y2 = sp.Symbol('Y2')

# constants
g = sp.Symbol('g')
m_1 = sp.Symbol('m_1')
m_2 = sp.Symbol('m_2')
l_1 = sp.Symbol('l_1')
l_2 = sp.Symbol('l_2')

# Kinetic energy
E_k1 = m_1 * l_1 * l_1 * x1 * x1
E_k2 = m_2 * (x1 * x1 * l_1 * l_1 + y1 * y1 * l_2 * l_2 + 2 * x1 * y1 * l_2 * l_2 * sp.cos(y0 - x0))

# Potential energy
E_p1 = m_1 * g * l_1 * (1-sp.cos(x0))
E_p2 = m_2 * g * (l_1 * (1 - sp.cos(x0)) + l_2 * (1 - sp.cos(y0)))

# Langangian
T = E_k1 + E_k2
V = E_p1 + E_p2
L = T - V

# with respect to theta
a = (L.diff(x1).diff(t) - L.diff(x0))
a = a.subs(x0.diff(t),X1).subs(y0.diff(t),Y1)
a = a.subs(X1.diff(t),X2).subs(Y1.diff(t),Y2)

# with respect to phi
b = (L.diff(y1).diff(t) - L.diff(y0))
b = b.subs(x0.diff(t),X1).subs(y0.diff(t),Y1)
b = b.subs(X1.diff(t),X2).subs(Y1.diff(t),Y2)

# solve and covert to list
eqa = sp.Eq(a,0)    # a = 0
eqb = sp.Eq(b, 0)   # b = 0
sol_dict = sp.solve([eqa, eqb], X2, Y2) # solve for x'' and y''
sol_keys = list(sol_dict.keys())
sol_vals = list(sol_dict.values())

print(sol_keys)
print(sol_vals)
