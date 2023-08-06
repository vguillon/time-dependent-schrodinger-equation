import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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

# Update scheme for FuncAnimation
def update (frame, normalize, schrodinger, x, x0_particle, title, line1_1, line2_1, line3_1):
    _, psi0 = schrodinger.psi0(x0_particle, particle['k0'], particle['sigma0'])
    psi_t = schrodinger.update(psi0, frame*simulation['T']/normalize) # Normalise should be equal to frame
    # |psi|^2
    psi_t2 = psi_t*np.conj(psi_t)
    psi_t2 = psi_t2.real
    line1_1.set_data(x, psi_t2)
    # Re(psi)
    re_psi_t = psi_t.real
    line2_1.set_data(x, re_psi_t)
    # Im(psi)
    im_psi_t = psi_t.imag
    line3_1.set_data(x, im_psi_t)
    # Display time
    title.set_text(f"t = {frame/normalize} T")
    return line1_1, line2_1, line3_1

# Setup psi, potential and figure. Then animate with the function update above
def animate (N, format, x0_particle):
    N = int(N)
    format = format.lower()
    x0_particle = float(x0_particle)
    # Intanciation
    schrodinger = Schrodinger(N, simulation['L'])
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
    # Plot setup
    fig = plt.figure(figsize=(8,10))
    ax1, ax2, ax3 = fig.subplots(3)
    x = np.linspace(0.0, simulation['L'], N)
    # For |psi|^2
    ax1.set_xlim((0.0, 1.0))
    ax1.set_ylim((-0.1, 3.0))
    line1_1, = ax1.plot([], [], lw=1, color='b', label=r'$\vert\psi(x,t)\vert^{2}$')
    line1_2, = ax1.plot([], [], lw=1, color='k', label=r'$V(x)$')
    line1_2.set_data(x, V)
    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$\vert\psi(x,t)\vert^{2}$')
    title = ax1.set_title(r'') # Display time (in update function)
    ax1.legend()
    # For real part of psi
    ax2.set_xlim((0.0, 1.0))
    ax2.set_ylim((-2.0, 2.0))
    line2_1, = ax2.plot([], [], lw=1, color='r', label=r'$\text{Re}~\psi(x,t)$')
    line2_2, = ax2.plot([], [], lw=1, color='k', label=r'$V(x)$')
    line2_2.set_data(x, V)
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$\text{Re}~\psi(x,t)$')
    ax2.legend()
    # For imaginary part of psi
    ax3.set_xlim((0.0, 1.0))
    ax3.set_ylim((-2.0, 2.0))
    line3_1, = ax3.plot([], [], lw=1, color='g', label=r'$\text{Im}~\psi(x,t)$')
    line3_2, = ax3.plot([], [], lw=1, color='k', label=r'$V(x)$')
    line3_2.set_data(x, V)
    ax3.set_xlabel(r'$x$')
    ax3.set_ylabel(r'$\text{Im}~\psi(x,t)$')
    ax3.legend()
    # Animation
    frames_number = 100
    ani = FuncAnimation(fig, update, frames=frames_number,
                        fargs=(frames_number, schrodinger, x, x0_particle,
                               title, line1_1, line2_1, line3_1,),
                               interval=50, blit=True)
    ani.save(f"{potential}_potential.{format}", fps=30)


