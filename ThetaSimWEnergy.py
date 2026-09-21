import math as m
import numpy as np
import matplotlib.pyplot as plt
theta_0 = 0.1
omega_0 = 0
g_da = 9.82
length = 1
mass = 1

h = 0.001
t_end = 10

# angular acceleration
def alphafunc(y,g,l):
    return -m.sin(y)*g/l

def StDevPop(list_obj): # Standard deviation for a population, not a sample
    mean = np.mean(list_obj)
    stdev = np.std(list_obj,ddof=0)
    return mean,stdev

def PlotStDev(list_obj,sigmawidth=5,res=100):
    mean,stdev = StDevPop(list_obj)
    var = stdev*stdev
    start = mean - sigmawidth*stdev
    end = mean + sigmawidth*stdev
    x_list = np.linspace(start,end,res)
    print(x_list)
    p = lambda x: m.exp(-((x-mean)**2/(2*var)))/(m.sqrt(2*m.pi)*stdev)
    y_list = list(map(p,x_list))
    plt.plot(x_list,y_list)
    plt.show()


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

def ThetaEnergyShiftWave(theta_list, Energy_list):
    theta_list = np.array(theta_list)/max(theta_list)
    Energy_list = np.array(Energy_list)/max(Energy_list)
    if len(theta_list) != len(Energy_list):
        raise ValueError('Theta and Energy list must have same length')
    wave = list(map(lambda a,b: a-b, theta_list, Energy_list))
    return wave




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

wave = ThetaEnergyShiftWave(theta_list, Energy_list)
plt.plot(t_list,wave)
plt.show()

#print(StDevPop(Energy_list))
#PlotStDev(Energy_list,5)