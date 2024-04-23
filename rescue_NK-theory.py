import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import math
import scipy.integrate as integrate
import scipy.special as sc
from scipy.optimize import root


# PARAMETERS
Wmax = 1.5			#maximum reproductive value
K  = 10000	                #carry capacity
N0 = 10000	                #initial pop size
mut = 0.001	                #mutation rate
L = 12                          # genome size
b = 6                          #range of f-values

C = np.sqrt(3*L/(np.pi*b*b))



############################################################################################################################################
# [Functions]
############################################################################################################################################


def p(u):
    return C*math.exp(-3*L*u*u/(b*b))

def fitness(u):
    return math.exp(-np.power((u-Delta),2)/2) * (Wmax) 

def prob_fix(u):
    pi = 1.0
    precisao = 1.0
    while precisao>1.e-14:
        pi_next = 1 - math.exp(-pi*fitness(u))
        precisao=math.fabs(pi_next-pi)
        pi = pi_next
    return (pi_next)



############################################################################################################################################
# [Integration]
############################################################################################################################################


# [ Integration in mutation (u) ]
def integrand_u(u): 
    return (p(u)*np.exp(-prob_fix(u)*A/L))		### genotypic
def integrand_frac_u(u): 
    return (p(u))



############################################################################################################################################
# [MAIN]
############################################################################################################################################

DEL = [0.01,0.02,0.03,0.05,0.07,0.08,0.1,0.12,0.15,0.17,0.2,0.22,0.25,0.27,0.3,0.32,0.35,0.37,0.4,0.42,0.45,0.47,0.5,0.52,0.55,0.57,0.6,0.62,0.65,0.67,0.70,0.72,0.75,0.77,0.80]
for delta in DEL:

 Delta = np.sqrt(-2*np.log((1-delta)/Wmax))
 A  = mut*N0*(1-delta)/delta
 umin = Delta - np.sqrt(2*np.log(Wmax))
 umax = Delta + np.sqrt(2*np.log(Wmax))
     
 frac = integrate.quad(integrand_frac_u, umin, umax)
 pi = integrate.quad(integrand_u, umin, umax)

 
 R = frac[0] - pi[0]
 Pr = (1 - R)**L
 
 f = open("theory_NK_gen-U%.5f.txt"%mut, "a")
 print("%.2f"%delta,"   ",Pr,"   ",pi[0],file=f)



f.close()



