import numpy as np 

def load_potential_parameters():
    potential, V0, width, x0, sigma = np.loadtxt('potential_parameters.txt', dtype=str, unpack=True)
    return potential, V0, width, x0, sigma

particle = {'k0': 2*np.pi/0.02,
             'sigma0': 0.04,
             'v': 4*np.pi/0.02, # v=2*k0
             'E': (2*np.pi/0.02)**2 # E = (k0)**2 if hbar=1 and m=1/2
             } 

simulation = {'L': 1.0, # lenght along x-axis
              'T': 1.0/particle['v'] # simulation 'time'
              }


