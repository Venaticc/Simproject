import math

from math import *

from sympy import root

Y_0= 10
L= 24
V_0= 0
h= 0.1
g= -9.82
a_0 = 0
O_0 = acos((L-Y_0)/(L))
B = abs(O_0)
x_0 = math.sin(O_0)*L
a_0 = (-sin(O_0)*g)/L
print(B)
omega = (a_0*h)/L
for i in range (1000000):
    if  abs(O_0 -((omega)*h))>O_0:
        a_0 = (-sin(O_0) * g) / L
        omega = (a_0 * h) / L
    else:
        O_0 = (O_0 -((omega)*h))
        omega = -(((2*g)/L)*(((cos((O_0))-cos((B)))))**(1/2))
    print(O_0,omega)
print(omega,O_0)
#*copysign(1, O_0)







