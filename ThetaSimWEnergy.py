import math as m
import numpy as np
import matplotlib.pyplot as plt
theta_0 = 0.1
omega_0 = 0
g_da = 9.82
length = 0.10
mass = 1

h = 0.001
t_end = 10

# angular acceleration
def alphafunc(y,g,l):
    return -m.sin(y)*g/l

def Energy_M(theta_list, omega_list, g_da, l, mass):
    Energy_list = list(map(lambda theta, w: mass*l*(0.5*l*w**2 + g_da*(1 - m.cos(theta))), theta_list, omega_list))
    return Energy_list

# runs the sim
def RunThetaSim(theta_0, t_end, h, omega, g_da, length):
    num_steps = int(t_end/h)                    # steps
    theta_list = list(np.zeros(num_steps+1))    # list of angles
    theta_list[0] = theta_0                     # saves first angle
    omega_list = list(np.zeros(num_steps+1))    # list of angular velocities
    omega_list[0] = omega_0                     # saves first angular velocity
    for i in np.arange(int(t_end/h)):           # recursive step siM
        theta = theta_list[i]                   # gets angle
        alpha = alphafunc(theta,g_da,length)    # current angular acceleration
        omega = omega + alpha*h                 # angle velocity
        theta = theta+omega*h                   # angle
        omega_list[i+1] = omega                 # saves angular velocity
        theta_list[i+1] = theta                 # saves angle

    t_list = np.linspace(0,t_end,num_steps+1)   # time_stamps

    return theta_list, omega_list, t_list

theta_list, omega_list, t_list = RunThetaSim(theta_0, t_end, h, omega_0, g_da, length)
Energy_list = Energy_M(theta_list, omega_list, g_da, length, mass)
print(theta_list)
print(Energy_list)
plt.subplot(2,1,1)
plt.plot(t_list,theta_list)
plt.title('Theta Simulation')
plt.xlabel('t')
plt.ylabel('theta')
plt.subplot(2,1,2)
plt.plot(t_list,Energy_list)
plt.title('Energy Simulation')
plt.xlabel('t')
plt.ylabel('energy')
plt.show()