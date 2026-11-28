import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# constants
g = sp.Symbol('g')
m_1 = sp.Symbol('m_1')
m_2 = sp.Symbol('m_2')
l_1 = sp.Symbol('l_1')
l_2 = sp.Symbol('l_2')

# assign values to constants
symbolvaluedict = {g: 9.82, m_1: 10, m_2: 1, l_1: 3, l_2: 1}

def EnergyErrorSTD(energy_array, sigma_width = 5, res = 100):
    mean = np.mean(energy_array)
    stdev = np.std(energy_array,ddof=0)
    if sigma_width == 0:
        return mean, stdev
    resolution = np.linspace(mean - sigma_width * stdev, mean + sigma_width * stdev, res)
    pdf = lambda x: m.exp(-((x-mean)**2/(2*stdev*stdev)))/(m.sqrt(2*m.pi)*stdev)
    pdf_map = list(map(pdf,resolution))
    return resolution, pdf_map, mean, stdev

def EnergyDeviation(energy_array):
    zeroed_array = energy_array - energy_array[0]
    zeroed_array *= zeroed_array
    return zeroed_array, np.sum(zeroed_array)

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
    # E_p3 = m_3 * g * (l_1 * (1 - sp.cos(x0)) + l_2 * (1 - sp.cos(y0)) + l_3 * (1 - sp.cos(z0)))

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

E_m, X0, Y0, X1, Y1, X2, Y2 = DeriveFunctions(symbolvaluedict)


X2lambda = sp.lambdify((X0, Y0, X1, Y1), X2)

Y2lambda = sp.lambdify((X0, Y0, X1, Y1), Y2)

E_mlambda = sp.lambdify((X0, Y0, X1, Y1), E_m)

def CalculateESM(h, t, angleX, angleY, omegaX, omegaY):
    omegaX = h*X2lambda(angleX, angleY, omegaX, omegaY)
    omegaY = h*Y2lambda(angleX, angleY, omegaX, omegaY)

    angleX = omegaX*h
    angleY = omegaY*h

    return angleX, angleY, omegaX, omegaY


def SimulateDPESM(h, t_start, t_end, angleX, angleY, omegaX, omegaY):
    num_steps = int((t_end - t_start)/h)
    t = t_start

    states = np.zeros((num_steps + 1, 4))
    states[0] = np.array([angleX, angleY, omegaX, omegaY])

    energy_array = np.zeros(num_steps + 1)
    energy_array[0] = E_mlambda(angleX, angleY, omegaX, omegaY)

    for i in np.arange(num_steps):
        states[i + 1] = CalculateESM(h, t, *states[i])

        t = t + h

        energy_array[i + 1] = E_mlambda(*states[i + 1])

    t_array = np.linspace(t_start, t_end, num_steps + 1, endpoint=True)

    return t_array, energy_array, *states.T

t_array, E_m_array, anX_sim, anY_sim, omX_sim, omY_sim = SimulateDPESM(0.001, 0, 20, m.pi, 2, 0, 0)

plt.plot(t_array, anX_sim, '-b')
plt.plot(t_array, anY_sim, '-r')
plt.show()

res, pdf_array, mean, stdev = EnergyErrorSTD(E_m_array)
print(mean, stdev)

plt.plot(res, pdf_array)
plt.show()

Edsq_array, E_dev = EnergyDeviation(E_m_array)
print(E_dev)

plt.plot(t_array, Edsq_array)
plt.show()

plt.plot(t_array, E_m_array, 'g-')
plt.show()
print(E_m_array)
print(omX_sim)
print(omY_sim)