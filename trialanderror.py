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
symbolvaluedict = {g: 9.82, m_1: 20, m_2: 2, l_1: 1, l_2: 3}

t = sp.Symbol('t')

X0 = sp.Function('X0')(t)
X1 = sp.Function('X1')(t)  # replacement function, magic trick
X2 = sp.Symbol('X2')

Y0 = sp.Function('Y0')(t)
Y1 = sp.Function('Y1')(t)  # also magic tricks
Y2 = sp.Symbol('Y2')

Eq1 = (m_1 + m_2)*l_1*l_1*X2 + m_2*l_1*l_2*Y2*sp.cos(Y0 - X0) - m_2 * l_1 * l_2 * Y1 * Y1 * sp.sin(Y0 - X0) + (m_1 + m_2) * g * l_1 * sp.sin(X0)

Eq2 = m_2 * l_2 * l_2 * Y2 + m_2 * l_1 * l_2 * X2 * sp.cos(Y0 - X0) + m_2 * l_1 * l_2 * X1 * X1 * sp.sin(Y0 - X0) + m_2 * g * l_2 * sp.sin(Y0)

sol_dict = sp.solve([Eq1, Eq2], X2, Y2)

X2 = sp.simplify(sol_dict[X2]).subs(symbolvaluedict)
Y2 = sp.simplify(sol_dict[Y2]).subs(symbolvaluedict)

print(X2.expand())
print(Y2.expand())

E_m = g*l_1*m_1*(1 - sp.cos(X0)) + g*m_2*(l_1*(1 - sp.cos(X0)) + l_2*(1 - sp.cos(Y0))) + l_1**2*m_1*X1**2/2 + m_2*(l_1**2*X1**2 + 2*l_2**2*X1*Y1*sp.cos(X0 - Y0) + l_2**2*Y1**2)/2
E_m = E_m.subs(symbolvaluedict)


def X2func(anX, anY, omX, omY):
    return X2.subs({X0:anX, Y0:anY, X1:omX, Y1:omY}).evalf()

def Y2func(anX, anY, omX, omY):
    return Y2.subs({X0:anX, Y0:anY, X1:omX, Y1:omY}).evalf()

def E_mfunc(anX, anY, omX, omY):
    return E_m.subs({X0:anX, Y0:anY, X1:omX, Y1:omY}).evalf()

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

        alphaX = X2func(anX, omX, anY, omY)
        alphaY = Y2func(anX, omX, anY, omY)

        omX = omX + alphaX * h
        omY = omY + alphaY * h

        anX = anX + omX * h
        anY = anY + omY * h

        anX_array[i+1] = anX % (2*m.pi)
        anY_array[i+1] = anY % (2*m.pi)
        omX_array[i+1] = omX
        omY_array[i+1] = omY

        energy_array[i+1] = E_mfunc(anX, omX, anY, omY)

    t_array = np.linspace(0, t_end, num_steps + 1)

    return t_array, anX_array, anY_array, omX_array, omY_array, energy_array

t_array, anX_sim, anY_sim, omX_sim, omY_sim, energy_sim = RunDoublePendulum(1, 1, 0, 0, 0.01, 20)

plt.plot(t_array, anX_sim, 'b-', label='P1')
plt.plot(t_array, anY_sim, 'r-', label='P2')
plt.show()

plt.plot(t_array, energy_sim, 'g-')
plt.show()


print(anX_sim)
print(anY_sim)
print(t_array)