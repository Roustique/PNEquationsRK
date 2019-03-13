# PNEquationsRK

The Runge-Kutta-4 Integrator for Post-Newtonian and Newtonian Equations of Motion

## Installing

Copy all files to your computer, then run
```
make
```
One should also create a directories ./input and ./output
## Instructions

### Data representation

Program reads data from file "./input/data.dat" and writes it to "./output/result.dat".
Input looks like:
```
Time step value (dt) (in c/(10 km) units)
Number of iterations (n)
first component of initial radius-vector (x) (in units of 10 km)
second component of initial radius-vector (y) (in units of 10 km)
first component of initial vector of speed (vx) (in units of c)
first component of initial vector of speed (vy) (in units of c)
mass of central body (M) (in units of Solar mass)
```

### Usage

Run compiled program with key "-P" to use Post-Newtonian Approximation
```
./main.out -P
```
Run it with "-N" to use Newtonian Approximation
```
./main.out -N
```
A Python Script named PlotDrawer.p could be used to create plots and also to compute some parameters, such as periapse precession, approximated semi-major axis and eccentricity for post-newtonian trajectories

MultyRun.py provides running program multiple times with different initial values, with no need of rewriting them each time. To use it, the different arrays of initial values (dt,n,x,y,vx,vy,M) should be written in string with delimiter ", " in file "./input/dataMulty.dat". The script copies this arrays from each string to the file ./input/data.dat and then runs program and PlotDrawer script.
## Authors

* **Rustam Gainutdinov** - [Roustique](https://github.com/Roustique)

