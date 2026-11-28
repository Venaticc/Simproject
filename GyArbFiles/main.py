from tokenize import triple_quoted

import JacobianAttemptDoublePendulumRK4 as JADPRK4
import LyapunovExponentSymbolSim as LESS
import RungeKuttaDoublePendulum as DPRk4
import SympySupportedPendulum as DPESM
import SympySupportedSinglePendulum as SPESM

import numpy as np
import sympy

h = 0.001
t_start = 0
t_end = 20

t_types = np.array([])

t_q = 5
epsilon = 1e-8
renorm_every = 10

angleX = 2
angleY = -1
omegaX = 0
omegaY = 0

states = np.array([angleX, angleY, omegaX, omegaY])

def Sim(type, svd, t_array, states, t_q, epsilon, renorm_every):
    match type:
        case "":
            pass

