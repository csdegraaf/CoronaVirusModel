import numpy as np
"""

dS / dt = -beta * S * I / N
dE / dt = +beta * S * I / N - sigma * E
dI / dt = +sigma * E - gamma * I + c * R * I / N
dR / dt = gamma * I - c * R * I / N



yprime = [dS / dt  dE / dt dI / dt   dRdt]
'

input:
  t current time
  y vector of current soln values
  y(1) = S, y(2) = E, y(3) = I, y(4) = R

parameters in "params"
  beta, N, sigma, gamma, c, R_zero_array(table of values)

output: (col vector)
  yprime(1) = dS / dt
  yprime(2) = dE / dt
  yprime(3) = dI / dt
  yprime(4) = dR / dt

"""


def seir_function(t, y, params):
    
    R_zero_array = params.r_zero
    
    min_t = np.min(R_zero_array[:, 0])
    n_table = np.shape(R_zero_array)
    max_t = np.max(R_zero_array[:, 0])
    t_val = max(min_t, min(t, max_t))
    
    R_zero = np.interp(t_val, R_zero_array[:, 0], R_zero_array[:, 1])
    
    gamma = params.gamma
    
    beta = R_zero * gamma
    
    N = params.N
    sigma = params.sigma
    c = params.c
    
    S = y[0]
    E = y[1]
    I = y[2]
    R = y[3]
    
    yprime = np.zeros(4)
    
    yprime[0] = -beta * S * I / N
    yprime[1] = +beta * S * I / N - sigma * E
    yprime[2] = +sigma * E - gamma * I + c * R * I / N
    yprime[3] = gamma * I - c * R * I / N
    return yprime