# SEIR Model for Corona
This project can be used to model virus infections throughout a country. The code is an exact copy of the matlab code 
from the [website of Peter Forsyth](https://cs.uwaterloo.ca/~paforsyt/SEIR.html). The parameters that are being modeled 
are 

* S = susceptible population
* E = Exposed (infected, not yet infectious)
* I = Infectious (now can infect others)
* R = Removed (got sick, now recovered and immune, or died :( )
* N = total population = (S + E + I + R)


note: added cRI/N term:  disease mutates and can cause reinfection or immunity lost. This assumes that mutated form 
jumps to Infected population. It can also be that a mutated form jumps to the Exposed population. For now, we assume 
c=0 (no mutation has been observed).
  
* dS/dt = -beta* S * I/N
* dE/dt = +beta* S* I/N - sigma* E
* dI/dt = +sigma* E -gamma* I + c* R* I/N
* dR/dt = gamma* I -c* R* I/N
