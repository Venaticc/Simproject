import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

theta_0 = 0.5
phi_0 = -0.9
omega_0 = 0
rho_0 = 0
g_da = 9.82
length_1 = 0.3
length_2 = 0.1
mass_1 = 1
mass_2 = 2

h = 0.001
t_end = 5


def AlphaTheta(theta, constant):
    return m.sin(theta)*constant

def AlphaPhi(phi,constant):
    return m.sin(phi)*constant


def RunSim(theta_0, phi_0, omega_0, rho_0, g_da, length_1, length_2,mass_1, mass_2):
    num_steps = int(t_end/h)                    # steps

    theta_list = np.zeros(num_steps+1)          # list of angles
    theta_list[0] = theta_0                     # saves first angle
    omega_list = np.zeros(num_steps+1)          # list of angular velocities
    omega_list[0] = omega_0                     # saves first angular velocity

    phi_list = np.zeros(num_steps+1)
    phi_list[0] = phi_0
    rho_list = np.zeros(num_steps+1)
    rho_list[0] = rho_0

    constant_alphafunc = -g_da*length_1*(mass_1 + mass_2)/(mass_1*length_1*length_1 + mass_2*length_2*length_2)
    constant_phifunc = -g_da/length_2

    print(constant_alphafunc)
    print(constant_phifunc)

    for i in np.arange(int(t_end/h)):           # recursive step siM

        theta = theta_list[i]                   # gets angle
        omega = omega_list[i]
        aTheta = AlphaTheta(theta,constant_alphafunc)    # current angular acceleration
        omega = omega + aTheta*h                 # angle velocity
        theta = theta+omega*h                   # angle

        omega_list[i+1] = omega                 # saves angular velocity
        theta_list[i+1] = theta                 # saves angle

        phi = phi_list[i]
        rho = rho_list[i]
        aPhi = AlphaPhi(phi,constant_phifunc)
        rho = rho + aPhi*h
        phi = phi + rho*h

        rho_list[i+1] = rho
        phi_list[i+1] = phi

    t_list = np.linspace(0,t_end,num_steps+1)   # time_stamps

    return theta_list, phi_list, (omega_list, rho_list), t_list


theta_list, phi_list, omega_rho_tuple, t_list = RunSim(theta_0, phi_0, omega_0, rho_0, g_da, length_1, length_2, mass_1, mass_2)

plt.subplot(2,1,2)
plt.plot(t_list, theta_list, 'b-', label='theta')
plt.subplot(2,1,1)
plt.plot(t_list, phi_list, 'r-', label='phi')
plt.show()

P_1 = list(map(lambda theta:  length_2 + length_1*(1- m.cos(theta)), theta_list))
P_2 = list(map(lambda P1y, phi: length_1 + length_2 - (P1y + length_2*m.cos(phi)), P_1, phi_list))
plt.plot(t_list, P_1, 'b-', label='P1')
plt.plot(t_list, P_2, 'r-', label='P2')
plt.show()