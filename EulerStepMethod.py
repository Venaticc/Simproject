import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math as m

# from sympy.abc import x, y

def carc(q):
    return 1/m.sqrt(q**2 + 1)

length = 1
g = 9.82

theta_0 = 0.97
b = carc(theta_0)

def fx(x,y):
    print(x,y,b)
    return -m.sqrt(2*g/length*(carc(x/y)-b))*carc(x/y)
def fy(x,y):
    return -m.sqrt(2*g/length*(carc(x/y)-b))*carc(x/y)*x/y
def ftheta(theta):
    return -m.sqrt(2*g/length*(carc(theta)-b))

x = length*b
y = length*b*theta_0
t_end = m.sqrt(2)+0.01
delta_t = 0.001

def Runftheta(theta):
    theta = theta - 0.0001
    theta_list = []
    for i in range(int(t_end/delta_t)):
        print(i)
        print(m.cos(theta), b)
        theta = theta + delta_t * ftheta(theta)
        theta_list.append(theta)
    x_list = list(map(lambda theta: length*m.cos(theta), theta_list))
    y_list = list(map(lambda theta: length*m.sin(theta), theta_list))
    print(x_list, y_list)
    plt.scatter(x_list, y_list)
    plt.show()

Runftheta(theta_0)






def Runfxy(x,y):
    x_list = []
    y_list = []
    for i in range(int(t_end/delta_t)):
        x = x + fx(x,y)*delta_t
        y = y + fy(x,y)*delta_t
        x_list.append(x)
        y_list.append(y)

    plt.scatter(x_list, y_list)