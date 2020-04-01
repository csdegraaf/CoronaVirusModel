"""
This code is a python version of the oringal Matlab/Octave code from Peter Forsyth
(see: https://cs.uwaterloo.ca/~paforsyt/SEIR.html)

parameters.py for SEIR model

S = susceptible population
E = Exposed (infected, not yet infectious)
I = Infectious (now can infect others)
R = Removed (got sick, now recovered and immune, or died :( )
N = total population = (S + E + I + R)


note: added cRI/N term:  disease
mutates, can cause reinfection, or immunity lost
This assumes that mutated form jumps to Infected population
Can also assume that mutated form jumps to Exposed population
For now, we assume c=0 (no mutation has been observed)

dS/dt = -beta*S*I/N
dE/dt = +beta*S*I/N - sigma*E
dI/dt = +sigma*E -gamma*I + c*R*I/N
dR/dt = gamma*I -c*R*I/N

this file passes seir_function in the calculations module to the ode solver
ode systen is specified in the calculations module

Parameters from:

Wang, H., Wang, Z., Dong, Y. et al. Phase-adjusted estimation of
the number of Coronavirus Disease 2019 cases in Wuhan, China.
Cell Discov 6, 10 (2020). https://doi.org/10.1038/s41421-020-0148-0

Wuhan, Feb 2020
Based on estimates for original outbreak in Wuhan

These p #rs are pretty much guestimates, but are probably
the right order of magnitude
"""
import numpy as np
import parameters as parameters
from scipy import integrate
from calculations_module import seir_function
import matplotlib.pyplot as plt


S_0 = 11.0e+6  # Wuhan city  excluding initial infected, exposed population,

I_0 = 40.0  # initial infected from market
# https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/newsâ€“wuhan-coronavirus/

E_0 = 20. * I_0  # initial exposed
# https://www.medrxiv.org/content/10.1101/2020.01.23.20018549v1.full.pdf


R_0 = 0  # initial recovered (not to be confused with R_zero, below)
# initially, no one has recovered

c = 0.0  # no mutation (yet)
# maybe this happens later?

N = S_0 + I_0 + E_0 + R_0  # N = total population

sigma = 1. / 5.2  # https://doi.org/10.1056/NEJMoa2001316 (2020).

gamma = 1. / 18.  # https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/newsâ€“wuhan-coronavirus
"""
 R_zero = number of people infected by each infectious person
          this has nothing to do with "R" = removed above
          or R_0 (initial value of recovered)
          but is common terminology (confusing, but usual notation)
     time dependent, starts offf large, than drops with
         time due to public health actions (i.e. quarantine, social distancing)

    R_zero > 1, cases increase
    R_zero < 1 cases peak and then drop off 

      R_zero declining with time https://www.nature.com/articles/s41421-020-0148-0

      beta = R_zero*gammma (done in "seir.m" )
 
     table of:   time(days)  R_zero
                  ....     ....
                  ....     ....
                  ....     ....
       linearly interpolate between times

       Note: this is different from Wang et al (2020), which assumes
             piecewise constant values for R_zero

"""
r_zero_array = np.zeros([6, 2])
r_zero_array[0, :] = [0.0,  3.0]# t=0 days    R_zero = 3.0
r_zero_array[1, :] = [20.0,  2.6]# t = 60 days R_zero = 2.6
r_zero_array[2, :] = [70.0,  1.9]# t = 70 days R_zero = 1.9
r_zero_array[3, :] = [84.0,  1.0]# t = 84 days R_zero = 1.0
r_zero_array[4, :] = [90.0,  .50]# t = 90 days R_zero = .50
r_zero_array[5, :] = [1000, .50]# t = 1000 days R_zero =.50

params = parameters.Params(c, N, sigma, gamma, r_zero_array)

t_0 = 0
tspan = np.linspace(t_0, 181, 180)  # time in days

y_init = np.zeros(4)
y_init[0] = S_0
y_init[1] = E_0
y_init[2] = I_0
y_init[3] = R_0


def seir_with_params(t, y):
    return seir_function(t, y, params)


r = integrate.ode(seir_with_params).set_integrator("dopri5")
r.set_initial_value(y_init, t_0)
y = np.zeros((len(tspan), len(y_init)))
y[0, :] = y_init  # array for solution
for i in range(1, 180):
    y[i, :] = r.integrate(tspan[i])
    if not r.successful():
        raise RuntimeError("Could not integrate")


fig, axes = plt.subplots(ncols=2)
axes[0].plot(tspan, y[:, 0], color="b", label="S: susceptible")
axes[1].plot(tspan, y[:, 1], color="r", label="E: exposed")
axes[0].set(xlabel="time (days)", ylabel="S: susceptible")
axes[1].set(xlabel="time (days)", ylabel="E: exposed")

axes[0].legend()
axes[1].legend()
plt.show()

fig, axes = plt.subplots(ncols=2)
axes[0].plot(tspan, y[:, 2], color="b", label="I: infectious")
axes[1].plot(tspan, y[:, 3], color="r", label="R: recovered")
axes[0].set(xlabel="time (days)", ylabel="I: infectious")
axes[1].set(xlabel="time (days)", ylabel="R: recovered")

axes[0].legend()
axes[1].legend()
plt.show()

total_cases = y[:, 1] + y[:, 2] + y[:, 3]
total_cases_active = y[:, 1] + y[:, 2]

fig, axes = plt.subplots(ncols=2)
axes[0].plot(tspan, total_cases, color="b", label="E+I+R: Total cases")
axes[1].plot(tspan, total_cases_active, color="r", label="E+I: Total cases")
axes[0].set(xlabel="time (days)", ylabel="E+I+R: Total cases")
axes[1].set(xlabel="time (days)", ylabel="E+I: Active cases")
plt.show()

nsteps = np.size(tspan)
S_end = y[nsteps - 1, 0]
E_end = y[nsteps - 1, 1]
I_end = y[nsteps - 1, 2]
R_end = y[nsteps - 1, 3]

total = S_end + E_end + I_end + R_end

print('\n time (days): ', tspan[nsteps-1], '\n')

print('total population: ', total)

print('initial infected: ', I_0)

print('\n total cases (E+I+R) at t= ', tspan[nsteps-1], ' \n', E_end + I_end + R_end)

print('\n Recovered at t= ', tspan[nsteps-1], ': ', R_end, ' \n')
print('\n Infected (infectious) at t= ', tspan[nsteps-1], ': ', I_end, ' \n')
print('\n Exposed (non-infectious) at t= ', tspan[nsteps-1], ': ', E_end, ' \n')
print('\n Susceptable at t= ', tspan[nsteps-1], ': ', S_end)
