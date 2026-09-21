import math as m
import numpy as np
import matplotlib.pyplot as plt
theta_0 = 0.1
omega_0 = 0
g_da = 9.82
length = 0.10

h = 0.001
t_end = 10

# angular acceleration
def alphafunc(y,g,l):
    return -m.sin(y)*g/l

# runs the sim
def RunThetaSim(theta_0, t_end, h, omega, g_da, length):
    num_steps = int(t_end/h)                    # steps
    theta_list = list(np.zeros(num_steps+1))    # list of angles
    theta_list[0] = theta_0                     # saves first angle
    for i in np.arange(int(t_end/h)):           # recursive step siM
        theta = theta_list[i]                   # gets angle
        alpha = alphafunc(theta,g_da,length)    # current angular acceleration
        omega = omega + alpha*h                 # angle velocity
        theta = theta+omega*h                   # angle
        theta_list[i+1] = theta                 # saves angle

    t_list = np.linspace(0,t_end,num_steps+1)   # time_stamps

    return theta_list, t_list

theta_list, t_list = RunThetaSim(theta_0, t_end, h, omega_0, g_da, length)
print(theta_list)
plt.plot(t_list,theta_list)
plt.show()