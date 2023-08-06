import numpy as np
import matplotlib.pyplot as plt 
from .crank_nicolson_algorithm import Schrodinger
from .parameters import load_potential_parameters, particle, simulation

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# Potential parameters
potential, V0, width, x0_potential, sigma = load_potential_parameters()
potential = str(potential).lower()
V0 = float(V0)
width = float(width)
x0_potential = float(x0_potential)
sigma = float(sigma)

def plot (N, format, t, x0_particle): # N, format, t and x0_particle are entered from the command line
    N = int(N)
    format = format.lower()
    t = float(t)
    x0_particle = float(x0_particle)
    # Intanciation
    schrodinger = Schrodinger(N, simulation['L'])
    # psi0
    x = np.linspace(0.0, simulation['L'], N)
    _, psi0 = schrodinger.psi0(x0_particle, particle['k0'], particle['sigma0'])
    psi02 = psi0*np.conj(psi0)
    psi02 = psi02.real # With the line above the imaganiray part is zero but still exists. We only keep the real part here so that matplotlib can make the plot
    # psi_t
    V = np.zeros(N)
    if potential == 'step':
        V = schrodinger.step_potential(V0)
    elif potential == 'barrier':
        V = schrodinger.barrier_potential(V0, width)
    elif potential == 'well':
        V = schrodinger.well_potential(V0, width)
    elif potential == 'gaussian':
        V = schrodinger.gaussian_potential(x0_potential, sigma)
    else:
        raise RuntimeError(f"{potential} is not a valid potential.")

    schrodinger.initialize(V*particle['E'])
    psi_t = schrodinger.update(psi0, t*simulation['T'])
    psi_t2 = psi_t*np.conj(psi_t)
    psi_t2 = psi_t2.real
    # Plots
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot()
    ax.plot(x, psi02, 'r', label=r'$|\psi(x, 0)|^{2}$')
    ax.plot(x, V, 'k', label=r'$V(x)$')
    ax.plot(x, psi_t2, 'b', label=r'$|\psi(x, t)|^{2}$')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$|\psi|^{2}$')
    ax.set_title(f"t = {t} T")
    ax.legend()
    fig.savefig(f"{potential}_potential.{format}", dpi=600, bbox_inches='tight')


