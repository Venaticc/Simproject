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
y=0
x=0
for i in range (1000000):
    if  (abs(O_0) -((omega)*h))>((B)) and O_0>0:
        omega = -omega
        print('x',y)
        y=+1
    elif (abs(O_0) +((omega)*h))>((B)) and O_0<0:
        omega = -omega
        print('y',x)
        x=x+1
    else:
        omega = omega
    O_0 = (O_0 -((omega)*h))
    omega = -(((2*g)/L)*(((cos((O_0))-cos((B)))))**(1/2))
    #print(O_0,omega)
    #print(abs(O_0) - ((omega) * h))
print(omega,O_0)
#*copysign(1, O_0)







