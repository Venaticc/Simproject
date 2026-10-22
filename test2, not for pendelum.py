import math as m
import numpy as np
import sympy as sp
from sympy.abc import x, y
import matplotlib.pyplot as plt

f = sp.Function('f')

g = sp.dsolve(f.diff(x,x) + 9*f(x),f(x))

print(g)