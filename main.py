from math import *
import matplotlib.pyplot as plt
x=[0]
y=[0.1]
w = 0
h = 0.001
g = 9.82
l = 0.10
O = ((0.1-(cos(y[0])*0.1))*g)
print(O)
W = [0]
Energi = [O]
def hogerlede(x,y,g,l):
    f = -sin(y)*g/l
    return f
B = 5/h
B = int(B)
for i in range(0,B):
    X = x[i]
    Y = y[i]
    dervitiv = hogerlede(X,Y,g,l)
    X = h+X
    w = w+dervitiv*h
    Y = Y+w*h
    W.append(w)
    if Y>0.1:
        print (Y,X)
    x.append(X)
    Energi.append(((0.1-cos(Y)*0.1)*g+((w*0.1)**2)/2))
    y.append(Y)
print(Y)
plt.plot(x,Energi)
plt.show()
plt.plot(x,W)
plt.show()
plt.plot(x,y)
plt.show()