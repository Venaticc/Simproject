import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import time

timestamp = time.time()

# constants
g = sp.Symbol('g')
m_1 = sp.Symbol('m_1')
m_2 = sp.Symbol('m_2')
l_1 = sp.Symbol('l_1')
l_2 = sp.Symbol('l_2')

# assign values to constants
svd = {g: 9.82, m_1: 10, m_2: 1, l_1: 3, l_2: 1}

# finds the functions for the angular acceleration of both angles and the total mechanical energy
t = sp.Symbol('t') # time

x0 = sp.Function('x0')(t)   # theta
x1 = x0.diff(t)             # theta prime
x2 = x1.diff(t)             # theta bis

X0 = sp.Function('X0')(t)
X1 = sp.Function('X1')(t)   # replacement function, magic trick
X2 = sp.Symbol('X2')

y0 = sp.Function('y0')(t)   # phi
y1 = y0.diff(t)             # phi prime
y2 = y1.diff(t)             # phi bis

Y0 = sp.Function('Y0')(t)
Y1 = sp.Function('Y1')(t)   # also magic tricks
Y2 = sp.Symbol('Y2')

# Kinetic energy E_k = mv^2/2
E_k1 = m_1 * l_1 * l_1 * x1 * x1 / 2
E_k2 = m_2 * (x1 * x1 * l_1 * l_1 + y1 * y1 * l_2 * l_2 + 2 * x1 * y1 * l_2 * l_2 * sp.cos(y0 - x0)) / 2

# Potential energy E_p = mgh
E_p1 = m_1 * g * l_1 * (1-sp.cos(x0))
E_p2 = m_2 * g * (l_1 * (1 - sp.cos(x0)) + l_2 * (1 - sp.cos(y0)))

from sympy import init_session
# Langangian
T = E_k1 + E_k2     # sum of E_k
print(T)
V = E_p1 + E_p2     # sum of E_p
print(V)
L = T - V           # the lagrangian
E_m = T + V         # sum of all energy

sp.init_session(quiet=True)