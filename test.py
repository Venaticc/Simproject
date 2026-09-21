import matplotlib.pyplot as plt
import numpy as np


def f(x,y):
    return 1+y*y

h = 0.001
y = 0
x = 0
x_end = 1.561

x_list = [x]
y_list = [y]
RK = False
for _ in range(int((x_end-x)/h)):
    if False:
        f_1 = f(x,y)
        f_2 = f(x+h/2,y+f_1*h/2)
        f_3 = f(x+h/2,y+f_2*h/2)
        f_4 = f(x+h,y+f_3*h)
        y = y+h*(f_1+2*f_2+2*f_3+f_4)/6
    if True:
        y = y + f(x,y)*h
    x = x+h
    x_list.append(x)
    y_list.append(y)

print(x_list[-2])
print(y_list[-2])
print(y_list)
print(x_list)
plt.plot(x_list, y_list)
plt.show()

"""
Lösa diffekvationer (först matematisk pendel och sedan dubbelpendel) med stegmetoderna Euler, Runge-Kutta och bakvänd Euler
Beräkna energis bevarande för att på så sätt bestämma hur verklighetstrogna de är
Potentiellt, om tid finns, pröva att bestämma Ljapunov-exponenten, som är ett mått på kaos


"""