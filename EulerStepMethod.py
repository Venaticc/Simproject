import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math as m

from sympy.abc import x, y

def carc(q):
    return 1/m.sqrt(q**2 + 1)

length = 1
g = 9.82
b = carc(12)

fx = m.sqrt(2*g/length*(carc(x/y)-b))*
fy = sym.diff(Fy, y)
print(fx, fy)

t_start = 0
t_end = 1
for  in range(t_start, t_stop)