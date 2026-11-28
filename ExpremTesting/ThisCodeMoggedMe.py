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
symbolvaluedict = {g: 9.82, m_1: 10, m_2: 2, l_1: 3, l_2: 4}


# the Lyapunov-exponent
def LyapunovExponent(omega_array):
    return np.mean(np.log(abs(omega_array)))


# finds the functions for the angular acceleration of both angles and the total mechanical energy
def DeriveFunctions(svd):
    t = sp.Symbol('t') # time

    x0 = sp.Function('x0')(t)   # theta
    x1 = x0.diff(t)             # theta prime

    X0 = sp.Function('X0')(t)
    X1 = sp.Function('X1')(t)   # replacement function, magic trick
    X2 = sp.Symbol('X2')

    y0 = sp.Function('y0')(t)   # phi
    y1 = y0.diff(t)             # phi prime

    Y0 = sp.Function('Y0')(t)
    Y1 = sp.Function('Y1')(t)   # also magic tricks
    Y2 = sp.Symbol('Y2')

    # Kinetic energy E_k = mv^2/2
    E_k1 = m_1 * l_1 * l_1 * x1 * x1 / 2
    E_k2 = m_2 * (x1 * x1 * l_1 * l_1 + y1 * y1 * l_2 * l_2 + 2 * x1 * y1 * l_2 * l_1 * sp.cos(y0 - x0)) / 2

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
    sol_dict = sp.solve([x_expr, y_expr], X2, Y2) # solve for x'' and y''

    # Finished functions!!!
    X2 = sp.simplify(sol_dict[X2].subs(svd))
    Y2 = sp.simplify(sol_dict[Y2].subs(svd))

    E_m = sp.simplify(E_m.subs(svd))

    return E_m, X0, Y0, X1, Y1, X2, Y2


# creates the functions
E_m, X0, Y0, X1, Y1, X2, Y2 = DeriveFunctions(symbolvaluedict)

E_mlambda = sp.lambdify((X0, Y0, X1, Y1), E_m, 'numpy')

accel = sp.lambdify((X0, Y0, X1, Y1), [X2, Y2], modules='math', cse=True)

def Derivative(a, b, c, d):
    x2, y2 = accel(a, b, c, d)
    return c, d, x2, y2

def CalculateRK4(h, t, a, b, c, d):
    k1 = Derivative(a, b, c, d)
    k2 = Derivative(a + h/2*k1[0], b + h/2*k1[1], c + h/2*k1[2], d + h/2*k1[3])
    k3 = Derivative(a + h/2*k2[0], b + h/2*k2[1], c + h/2*k2[2], d + h/2*k2[3])
    k4 = Derivative(a + h*k3[0],   b + h*k3[1],   c + h*k3[2],   d + h*k3[3])
    return tuple(s + h/6*(p + 2*q + 2*r + u)
                 for s, p, q, r, u in zip((a, b, c, d), k1, k2, k3, k4))

def SimDPRK4(h, t_start, t_end, angleX, angleY, omegaX, omegaY):
    num_steps = int((t_end - t_start)/h)
    t = t_start

    states = np.zeros((num_steps + 1, 4))
    states[0] = np.array([angleX, angleY, omegaX, omegaY])

    energy_array = np.zeros(num_steps + 1)
    energy_array[0] = E_mlambda(angleX, angleY, omegaX, omegaY)

    for i in np.arange(num_steps):
        states[i + 1] = CalculateRK4(h, t, *states[i])

        t = t + h

        energy_array[i + 1] = E_mlambda(*states[i + 1])

    t_array = np.linspace(t_start, t_end, num_steps + 1, endpoint=True)

    return t_array, energy_array, *states.T

print(time.time() - timestamp)
timestamp = time.time()

SimVals = np.array(SimDPRK4(0.001, 0, 20, -1, 2, 0, 0))

print(time.time() - timestamp)

t_array = SimVals[0]
anX_sim = SimVals[2]
anY_sim = SimVals[3]

plt.plot(t_array, anX_sim, '-b')
plt.plot(t_array, anY_sim, '-r')
plt.show()
