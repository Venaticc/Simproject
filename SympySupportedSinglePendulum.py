import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# constants
g = sp.Symbol('g')
m_1 = sp.Symbol('m_1')
l_1 = sp.Symbol('l_1')

# assign values to constants
symbolvaluedict = {g: 9.82, m_1: 10, l_1: 3}

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

    # Kinetic energy E_k = mv^2/2
    E_k1 = m_1 * l_1 * l_1 * x1 * x1 / 2

    # Potential energy E_p = mgh
    E_p1 = m_1 * g * l_1 * (1-sp.cos(x0))

    # Langangian
    T = E_k1     # sum of E_k
    V = E_p1     # sum of E_p
    L = T - V           # the lagrangian
    E_m = T + V         # sum of all energy

    # with respect to theta
    x_expr = (L.diff(x1).diff(t) - L.diff(x0))
    x_expr = x_expr.subs(x0.diff(t),X1)
    x_expr = x_expr.subs(X1.diff(t),X2)
    x_expr = x_expr.subs(x0,X0)

    # with respect to phi

    # replacing symbols
    E_m = E_m.subs({x0:X0, x1:X1})

    # solve and convert to list
    eqx = sp.Eq(x_expr,0)    # x_expr = 0
    sol_dict = sp.solve([eqx], X2) # solve for x'' and y''
    # sol_keys = list(sol_dict.keys())      # probably not nessecary
    sol_vals = list(sol_dict.values())

    # Finished functions!!!
    X2 = sp.simplify(sol_dict[X2].subs(svd))

    E_m = sp.simplify(E_m.subs(svd))

    return E_m, X0, X1, X2

E_m, X0, X1, X2 = DeriveFunctions(symbolvaluedict)

X2lambda = sp.lambdify((X0, X1), X2)

E_mlambda = sp.lambdify((X0, X1), E_m)

def RunDoublePendulum(angleX, angleY, omegaX, omegaY, h, t_end):
    num_steps = int(t_end/h)

    anX_array = np.zeros(num_steps + 1)
    anX_array[0] = angleX % (2*m.pi)
    anY_array = np.zeros(num_steps + 1)
    anY_array[0] = angleY % (2*m.pi)
    omX_array = np.zeros(num_steps + 1)
    omX_array[0] = omegaX
    omY_array = np.zeros(num_steps + 1)
    omY_array[0] = omegaY

    energy_array = np.zeros(num_steps + 1)
    energy_array[0] = E_mfunc(angleX, omegaX, angleY, omegaY)

    for i in np.arange(num_steps):

        anX = anX_array[i]
        anY = anY_array[i]
        omX = omX_array[i]
        omY = omY_array[i]

        alphaX = X2lambda(anX, omX, anY, omY)
        alphaY = Y2lambda(anX, omX, anY, omY)

        omX = omX + alphaX * h
        omY = omY + alphaY * h

        anX = anX + omX * h
        anY = anY + omY * h

        anX_array[i+1] = anX
        anY_array[i+1] = anY
        omX_array[i+1] = omX
        omY_array[i+1] = omY

        energy_array[i+1] = E_mlambda(anX, omX, anY, omY)

    t_array = np.linspace(0, t_end, num_steps + 1)

    return t_array, anX_array, anY_array, omX_array, omY_array, energy_array

def Position(anX_array, anY_array, svd):
    P1_lambda = lambda x: ((l_1*sp.cos(x)).subs(svd)).evalf(4)
    P2_lambda = lambda y: ((l_2*sp.cos(y)).subs(svd)).evalf(4)

    P1_array = np.array(list(map(P1_lambda, anX_array)))
    P2_array = np.array(list(map(P2_lambda, anY_array))) + P1_array

    return P1_array, P2_array

t_array, anX_sim, anY_sim, omX_sim, omY_sim, E_m_array = RunDoublePendulum(m.pi, 1, 0, 0, 0.01, 20)

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