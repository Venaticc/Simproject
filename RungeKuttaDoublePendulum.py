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
symbolvaluedict = {g: 9.82, m_1: 10, m_2: 1, l_1: 3, l_2: 1}

# finds the functions for the angular acceleration of both angles and the total mechanical energy
def DeriveFunctions(svd):
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

    # Langangian
    T = E_k1 + E_k2     # sum of E_k
    V = E_p1 + E_p2     # sum of E_p
    L = T - V           # the lagrangian
    E_m = T + V         # sum of all energy

    # with respect to theta
    x_expr = (L.diff(x1).diff(t) - L.diff(x0))
    x_expr = x_expr.subs(x0.diff(t),X1).subs(y0.diff(t),Y1)
    x_expr = x_expr.subs(X1.diff(t),X2).subs(Y1.diff(t),Y2)
    x_expr = x_expr.subs(x0,X0).subs(y0,Y0)

    # with respect to phi
    y_expr = (L.diff(y1).diff(t) - L.diff(y0))
    y_expr = y_expr.subs(x0.diff(t),X1).subs(y0.diff(t),Y1)
    y_expr = y_expr.subs(X1.diff(t),X2).subs(Y1.diff(t),Y2)
    y_expr = y_expr.subs(x0,X0).subs(y0,Y0)

    # replacing symbols
    E_m = E_m.subs({x0:X0, x1:X1, y0:Y0, y1:Y1})

    # solve and convert to list
    eqx = sp.Eq(x_expr,0)    # x_expr = 0
    eqy = sp.Eq(y_expr, 0)   # y_expr = 0
    sol_dict = sp.solve([eqx, eqy], X2, Y2) # solve for x'' and y''
    # sol_keys = list(sol_dict.keys())      # probably not nessecary
    sol_vals = list(sol_dict.values())

    # Finished functions!!!
    X2 = sp.simplify(sol_dict[X2].subs(svd))
    Y2 = sp.simplify(sol_dict[Y2].subs(svd))

    E_m = sp.simplify(E_m.subs(svd))

    return E_m, X0, Y0, X1, Y1, X2, Y2


# creates the functions
E_m, X0, Y0, X1, Y1, X2, Y2 = DeriveFunctions(symbolvaluedict)


# X_n+1 = X_n + h*fX
def fX(t, angleX, angleY, omegaX, omegaY):
    return omegaX


# omegaX_n+1 = omegaX_n + h*gX
gXlambda = sp.lambdify((X0, Y0, X1, Y1), X2)
def gX(t, angleX, angleY, omegaX, omegaY):
    return gXlambda(angleX, angleY, omegaX, omegaY)


# Y_n+1 = Y_n + h*fY
def fY(t, angleX, angleY, omegaX, omegaY):
    return omegaY


# omegaY_n+1 = omegaY_n + h*gY
gYlambda = sp.lambdify((X0, Y0, X1, Y1), Y2)
def gY(t, angleX, angleY, omegaX, omegaY):
    return gYlambda(angleX, angleY, omegaX, omegaY)


# total mechanical energy
E_mlambda = sp.lambdify((X0, Y0, X1, Y1), E_m)


# Find RK4 for a double pendulum
def CalculateRK4(h, t, angleX, angleY, omegaX, omegaY):
    k_1X = h*fX(t, angleX, angleY, omegaX, omegaY) # X
    k_1Y = h*fY(t, angleX, angleY, omegaX, omegaY) # Y

    r_1X = h*gX(t, angleX, angleY, omegaX, omegaY) # OMEGA X = h * ALPHA X
    r_1Y = h*gY(t, angleX, angleY, omegaX, omegaY) # OMEGA Y = h * ALPHA Y

    k_2X = h*fX(t+h/2, angleX + k_1X/2, angleY + k_1Y/2, omegaX + r_1X/2, omegaY + r_1Y/2)
    k_2Y = h*fY(t+h/2, angleX + k_1X/2, angleY + k_1Y/2, omegaX + r_1X/2, omegaY + r_1Y/2)

    r_2X = h*gX(t+h/2, angleX + k_1X/2, angleY + k_1Y/2, omegaX + r_1X/2, omegaY + r_1Y/2)
    r_2Y = h*gY(t+h/2, angleX + k_1X/2, angleY + k_1Y/2, omegaX + r_1X/2, omegaY + r_1Y/2)

    k_3X = h*fX(t+h/2, angleX + k_2X/2, angleY + k_2Y/2, omegaX + r_2X/2, omegaY + r_2Y/2)
    k_3Y = h*fY(t+h/2, angleX + k_1X/2, angleY + k_1Y/2, omegaX + r_1X/2, omegaY + r_1Y/2)

    r_3X = h*gX(t+h/2, angleX + k_2X/2, angleY + k_2Y/2, omegaX + r_2X/2, omegaY + r_2Y/2)
    r_3Y = h*gY(t+h/2, angleX + k_2X/2, angleY + k_2Y/2, omegaX + r_2X/2, omegaY + r_2Y/2)

    k_4X = h*fX(t+h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)
    k_4Y = h*fY(t+h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)

    r_4X = h*gX(t+h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)
    r_4Y = h*gY(t+h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)

    angleX = angleX + k_1X / 6 + k_2X / 3 + k_3X / 3 + k_4X / 6
    angleY = angleY + k_1Y / 6 + k_2Y / 3 + k_3Y / 3 + k_4Y / 6

    omegaX = omegaX + r_1X / 6 + r_2X / 3 + r_3X / 3 + r_4X / 6
    omegaY = omegaY + r_1Y / 6 + r_2Y / 3 + r_3Y / 3 + r_4Y / 6

    return angleX, angleY, omegaX, omegaY


# Run RK4 step simulation for a double pendulum
def SimulateDoublePendulumRK4(h, t_start, t_end, angleX, angleY, omegaX, omegaY):
    num_steps = int((t_end - t_start)/h)

    t = t_start

    anX_array = np.zeros(num_steps + 1)
    anX_array[0] = angleX % (2*m.pi)
    anY_array = np.zeros(num_steps + 1)
    anY_array[0] = angleY % (2*m.pi)
    omX_array = np.zeros(num_steps + 1)
    omX_array[0] = omegaX
    omY_array = np.zeros(num_steps + 1)
    omY_array[0] = omegaY

    energy_array = np.zeros(num_steps + 1)
    energy_array[0] = E_mlambda(angleX, angleY, omegaX, omegaY)

    for i in np.arange(num_steps):
        anX = anX_array[i]
        anY = anY_array[i]
        omX = omX_array[i]
        omY = omY_array[i]

        anX, anY, omX, omY = CalculateRK4(h, t, anX, anY, omX, omY)

        t = t + h

        anX_array[i + 1] = anX
        anY_array[i + 1] = anY
        omX_array[i + 1] = omX
        omY_array[i + 1] = omY

        energy_array[i + 1] = E_mlambda(anX, anY, omX, omY)

    t_array = np.linspace(t_start, t_end, num_steps + 1, endpoint=True)

    return t_array, anX_array, anY_array, omX_array, omY_array, energy_array


# time to define all functions and expressions
print(time.time() - timestamp)
timestamp = time.time()


# Run the simulation
t_array, anX_sim, anY_sim, omX_sim, omY_sim, E_m_array = SimulateDoublePendulumRK4(0.001, 0, 20, m.pi, 1, 0, 0)


# Time to run sim
print(time.time() - timestamp)


# graph of angles over time
plt.plot(t_array, anX_sim, '-b')
plt.plot(t_array, anY_sim, '-r')
plt.show()


# evergy over time, tight
plt.plot(t_array, E_m_array, '-b')
plt.axis("tight")
plt.show()


# an assortment of data
print("Starting E_m:", E_m_array[0])
print("E max - min:", max(E_m_array) - min(E_m_array))
print(np.std(E_m_array, ddof=0))


# scaled and absoluted energy to visualize deviation from start
plt.plot(t_array, np.abs(E_m_array - E_m_array[0]), '-r')
plt.axis("tight")
plt.show()


