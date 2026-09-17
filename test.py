import matplotlib.pyplot as plt
import numpy as np


def f(x,y):
    return -0.5*x*y

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