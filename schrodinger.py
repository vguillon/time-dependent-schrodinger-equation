import argparse
from schrodinger.plots import plot
from schrodinger.animation import animate

parser = argparse.ArgumentParser(description='A program for solving the time-dependent one dimensional Schrodinger equation')

parser.add_argument('-N', dest='points_number', type=int, default='2000', help='Number of points along x-axis (default = 2000)')
parser.add_argument('-o', dest='output', type=str, default='snapshot', help='output type: snapshot or movie (default = snapshot)')
parser.add_argument('-t', dest='time', type=float, default='0.3',
                    help='time in unit of T: must be between 0.0 and 1.0 (default = 0.3)') # default time for snapshot is t = 0.3*T
parser.add_argument('-sf', dest='snapshot_format', type=str, default='pdf', help='snapshot format: pdf, png, ... (default = pdf)')
parser.add_argument('-mf', dest='movie_format', type=str, default='mp4',
                    help='movie format: see matplotlib.animation documentation for supported format (default = mp4)')
parser.add_argument('-x', dest='initial_pos', type=float, default='0.3',
                    help='wave packet initial position: between 0.0 and 1.0 (default = 0.3)') # By default the initial wave packet is centered at x=0.3

args = parser.parse_args() 
#print(args)

###
N = int(args.points_number)
if N < 2000:
    raise ValueError(f"N must be >= 2000 (here {N} < 2000)")

out = ['snapshot', 'movie']
output = str(args.output).lower()
if output not in out:
    raise RuntimeError(f"Invalid option '{output}' (try 'snapshot' or 'movie')")

t = float(args.time)
if t < 0.0 or t >= 1.0:
    raise ValueError(f"time must be in [0.0, 1.0[ (here time = {t})")

snapshot_format = str(args.snapshot_format).lower()

movie_format = str(args.movie_format).lower()

x0_particle = float(args.initial_pos)
if x0_particle < 0.0 or x0_particle >= 1.0:
    raise ValueError(f"Initial particle position must be in [0.0, 1.0[ (here x0 = {x0_particle})")

###
try:
    if output == 'snapshot':
        plot(N, snapshot_format, t, x0_particle)
    if output == 'movie':
        animate(N, movie_format, x0_particle)
except RuntimeError as re:
    print(type(re), ":", re)
except ValueError as ve:
    print(type(ve), ":", ve)


