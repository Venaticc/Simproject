import matplotlib.pyplot as plt
import numpy as np


def f(x,y):
    return -0.5*x*y

def g(x,y):
    return

h = 0.01
y = 0.1
x = -5
x_end = 5

x_list = []
y_list = []
for i in range(int((x_end-x)/h)):
    y = y+h*f(x,y)
    x = x+h
    x_list.append(x)
    y_list.append(y)

plt.plot(x_list, y_list)
plt.show()

"""
Lösa diffekvationer (först matematisk pendel och sedan dubbelpendel) med stegmetoderna Euler, Runge-Kutta och bakvänd Euler
Beräkna energis bevarande för att på så sätt bestämma hur verklighetstrogna de är
Potentiellt, om tid finns, pröva att bestämma Ljapunov-exponenten, som är ett mått på kaos


"""