# PNEquationsRK

The Runge-Kutta-4 Integrator for Post-Newtonian and Newtonian Equations of Motion

## Installing

Copy all files to your computer, then run
```
make
```
## Instructions

### Data representation

Program reads data from file "data.dat" and writes it to "result.dat".
Input looks like:
```
Time step value (dt)
Number of iterations (n)
first component of initial radius-vector (x)
second component of initial radius-vector (y)
first component of initial vector of speed (vx)
first component of initial vector of speed (vy)
gravitational parameter for central body (GM)
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
## Authors

* **Rustam Gainutdinov** - [Roustique](https://github.com/Roustique)

