# Examples

The initial wave packet is a (non-normalized) Gaussian function. It evolves over a period T along the x-axis. A potential is placed on its path.

![movie gaussian potential](https://github.com/vguillon/time-dependent-schrodinger-equation/blob/main/movies/gaussian_potential.gif)

![snapshot step potential at t=0.3T](https://github.com/vguillon/time-dependent-schrodinger-equation/blob/main/snapshots/step_potential.png)

# Installation

The program needs [numpy](https://numpy.org) and [matplotlib](https://matplotlib.org) to run. For movie generation, see the
[matplotlib.animation](https://matplotlib.org/stable/api/animation_api.html) documentation ([ffmpeg](https://www.ffmpeg.org) recommended).

You can clone the project in your current repository with:

```
$ git clone https://github.com/vguillon/time-dependent-schrodinger-equation.git
```

# Parameters of the potentials

You can choose the parameters of the potentials directly in the **potential_parameters.txt** file. The parameters are:

- `Potential`: select potential. Available potentials are 'step', 'barrier', 'well' and 'gaussian' ;
- `V0`: potential height. Only used for 'step', 'barrier' and 'well' ;
- `width`: potential width. Only used for 'barrier' and 'well'. Hints: for 'barrier', `width` = 0.03 and for 'well', `width` = 0.1 give good results ;
- `x0`: center of the gaussian potential ;
- `sigma`: width of the gaussian potential.

# Running the program

The program allows several command line options. Each of these options possess a default value. Hence, after the potential parameters are configured, you can just run 

```
$ python3 schrodinger.py
```

to obtain a snapshot similare to the one above.

To show the available command line options, type:

```
$ python3 schrodinger.py -h
```

For instance, if you want a .png snapshot of the wave function after half a period for a wave packet starting at x=0.1 at t=0, you can write:

```
$ python3 schrodinger.py -sf png -t 0.5 -x 0.1
```

If you want to create a .mp4 movie over a period, type:

```
$ python3 schrodinger.py -o movie
```

If you want to create a .gif over a period with 3000 points along the x-axis for a wave packet starting at x=0.2 at t=0, write:

```
$ python3 schrodinger.py -o movie -mf gif -N 3000 -x 0.2
```

The generation of the datas can take few minutes.
