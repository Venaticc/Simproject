import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

timestamp = time.time()

# constants
g = sp.Symbol('g')
m_1 = sp.Symbol('m_1')
m_2 = sp.Symbol('m_2')
l_1 = sp.Symbol('l_1')
l_2 = sp.Symbol('l_2')

# assign values to constants
symbolvaluedict = {g: 9.82, m_1: 10, m_2: 2, l_1: 3, l_2: 4}


l_sum = symbolvaluedict[l_1] + symbolvaluedict[l_2]


h = 0.001
skip_size = 50


def DeriveFunctions(svd):
    t = sp.Symbol('t')  # time

    x0 = sp.Function('x0')(t)  # theta
    x1 = x0.diff(t)  # theta prime
    x2 = x1.diff(t)  # theta bis

    X0 = sp.Function('X0')(t)
    X1 = sp.Function('X1')(t)  # replacement function, magic trick
    X2 = sp.Symbol('X2')

    y0 = sp.Function('y0')(t)  # phi
    y1 = y0.diff(t)  # phi prime
    y2 = y1.diff(t)  # phi bis

    Y0 = sp.Function('Y0')(t)
    Y1 = sp.Function('Y1')(t)  # also magic tricks
    Y2 = sp.Symbol('Y2')

    # Kinetic energy E_k = mv^2/2
    E_k1 = m_1 * l_1 * l_1 * x1 * x1 / 2
    E_k2 = m_2 * (x1 * x1 * l_1 * l_1 + y1 * y1 * l_2 * l_2 + 2 * x1 * y1 * l_2 * l_2 * sp.cos(y0 - x0)) / 2

    # Potential energy E_p = mgh
    E_p1 = m_1 * g * l_1 * (1 - sp.cos(x0))
    E_p2 = m_2 * g * (l_1 * (1 - sp.cos(x0)) + l_2 * (1 - sp.cos(y0)))

    # Langangian
    T = E_k1 + E_k2  # sum of E_k
    V = E_p1 + E_p2  # sum of E_p
    L = T - V  # the lagrangian
    E_m = T + V  # sum of all energy

    # with respect to theta
    x_expr = (L.diff(x1).diff(t) - L.diff(x0))
    x_expr = x_expr.subs(x0.diff(t), X1).subs(y0.diff(t), Y1)
    x_expr = x_expr.subs(X1.diff(t), X2).subs(Y1.diff(t), Y2)
    x_expr = x_expr.subs(x0, X0).subs(y0, Y0)

    # with respect to phi
    y_expr = (L.diff(y1).diff(t) - L.diff(y0))
    y_expr = y_expr.subs(x0.diff(t), X1).subs(y0.diff(t), Y1)
    y_expr = y_expr.subs(X1.diff(t), X2).subs(Y1.diff(t), Y2)
    y_expr = y_expr.subs(x0, X0).subs(y0, Y0)

    # replacing symbols
    E_m = E_m.subs({x0: X0, x1: X1, y0: Y0, y1: Y1})

    # solve and convert to list
    eqx = sp.Eq(x_expr, 0)  # x_expr = 0
    eqy = sp.Eq(y_expr, 0)  # y_expr = 0
    sol_dict = sp.solve([eqx, eqy], X2, Y2)  # solve for x'' and y''
    # sol_keys = list(sol_dict.keys())      # probably not nessecary
    sol_vals = list(sol_dict.values())

    # Finished functions!!!
    X2 = sp.simplify(sol_dict[X2].subs(svd))
    Y2 = sp.simplify(sol_dict[Y2].subs(svd))

    E_m = sp.simplify(E_m.subs(svd))

    return E_m, X0, Y0, X1, Y1, X2, Y2


E_m, X0, Y0, X1, Y1, X2, Y2 = DeriveFunctions(symbolvaluedict)


def fX(t, angleX, angleY, omegaX, omegaY):
    return omegaX


gXlambda = sp.lambdify((X0, Y0, X1, Y1), X2)


def gX(t, angleX, angleY, omegaX, omegaY):
    return gXlambda(angleX, angleY, omegaX, omegaY)


def fY(t, angleX, angleY, omegaX, omegaY):
    return omegaY


gYlambda = sp.lambdify((X0, Y0, X1, Y1), Y2)


def gY(t, angleX, angleY, omegaX, omegaY):
    return gYlambda(angleX, angleY, omegaX, omegaY)


E_mlambda = sp.lambdify((X0, Y0, X1, Y1), E_m)


def QuickCalculateRK4(h, t, angleX, angleY, omegaX, omegaY):
    k_1X = h * fX(t, angleX, angleY, omegaX, omegaY)  # OMEGA X
    k_1Y = h * fY(t, angleX, angleY, omegaX, omegaY)  # OMEGA Y

    r_1X = h * gX(t, angleX, angleY, omegaX, omegaY)  # ALPHA X
    r_1Y = h * gY(t, angleX, angleY, omegaX, omegaY)  # ALPHA Y

    k_2X = h * fX(t + h / 2, angleX + k_1X / 2, angleY + k_1Y / 2, omegaX + r_1X / 2, omegaY + r_1Y / 2)
    k_2Y = h * fY(t + h / 2, angleX + k_1X / 2, angleY + k_1Y / 2, omegaX + r_1X / 2, omegaY + r_1Y / 2)

    r_2X = h * gX(t + h / 2, angleX + k_1X / 2, angleY + k_1Y / 2, omegaX + r_1X / 2, omegaY + r_1Y / 2)
    r_2Y = h * gY(t + h / 2, angleX + k_1X / 2, angleY + k_1Y / 2, omegaX + r_1X / 2, omegaY + r_1Y / 2)

    k_3X = h * fX(t + h / 2, angleX + k_2X / 2, angleY + k_2Y / 2, omegaX + r_2X / 2, omegaY + r_2Y / 2)
    k_3Y = h * fY(t + h / 2, angleX + k_1X / 2, angleY + k_1Y / 2, omegaX + r_1X / 2, omegaY + r_1Y / 2)

    r_3X = h * gX(t + h / 2, angleX + k_2X / 2, angleY + k_2Y / 2, omegaX + r_2X / 2, omegaY + r_2Y / 2)
    r_3Y = h * gY(t + h / 2, angleX + k_2X / 2, angleY + k_2Y / 2, omegaX + r_2X / 2, omegaY + r_2Y / 2)

    k_4X = h * fX(t + h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)
    k_4Y = h * fY(t + h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)

    r_4X = h * gX(t + h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)
    r_4Y = h * gY(t + h, angleX + k_3X, angleY + k_3Y, omegaX + r_3X, omegaY + r_3Y)

    angleX = angleX + k_1X / 6 + k_2X / 3 + k_3X / 3 + k_4X / 6
    angleY = angleY + k_1Y / 6 + k_2Y / 3 + k_3Y / 3 + k_4Y / 6

    omegaX = omegaX + r_1X / 6 + r_2X / 3 + r_3X / 3 + r_4X / 6
    omegaY = omegaY + r_1Y / 6 + r_2Y / 3 + r_3Y / 3 + r_4Y / 6

    return angleX, angleY, omegaX, omegaY


def QuickDPRK4(h, t_start, t_end, angleX, angleY, omegaX, omegaY):
    # Skips saving angular velocities

    num_steps = int((t_end - t_start) / h)
    t = t_start

    anX_array = np.zeros(num_steps + 1)
    anX_array[0] = angleX
    anY_array = np.zeros(num_steps + 1)
    anY_array[0] = angleY

    for i in np.arange(num_steps):
        anX_array[i + 1], anY_array[i + 1], omegaX, omegaY = QuickCalculateRK4(h, t, anX_array[i], anY_array[i], omegaX, omegaY)
        t += h

    return np.linspace(t_start, t_end, num_steps + 1, endpoint=True), anX_array, anY_array


t_array, anX_sim, anY_sim = QuickDPRK4(h, 0, 50, m.pi/2, m.pi, 0, 0)


anim_h = h*skip_size

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(autoscale_on=False, xlim=(-l_sum, l_sum), ylim=(-l_sum, l_sum))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], '.-', lw=1, ms=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

F1 = list(anX_sim)[::skip_size]

F2 = ([], [])
l_1 = symbolvaluedict[l_1]
l_2 = symbolvaluedict[l_2]

x1_data = l_1 * np.sin(F1)

y1_data = -l_1 * np.cos(F1)

f3 = list(anY_sim)[::skip_size]
f4 = ([], [])
x2_data = l_2 * np.sin(f3) + x1_data
y2_data = -l_2 * np.cos(f3) + y1_data

t_anim = t_array[::skip_size]

def animate(i):
    thisx = [0, x1_data[i], x2_data[i]]
    thisy = [0, y1_data[i], y2_data[i]]

    history_x = x2_data[:i+1]
    history_y = y2_data[:i+1]
    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)
    time_text.set_text(time_template % (t_anim[i]))
    return line, trace, time_text


ani = animation.FuncAnimation(
    fig, animate, len(t_anim), interval=anim_h * len(t_anim), blit=True)
ani.save("animation.gif", writer='pillow', fps=round(1/anim_h), dpi=100)
plt.show()
