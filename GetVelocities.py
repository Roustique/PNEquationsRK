#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import fsolve

inputfile = open("./input/initvel.dat")
a = float(inputfile.readline())
e = float(inputfile.readline())
M = float(inputfile.readline())
G = 1.0
c = 1.0
phi = G*M/(a*(1-e**2))

def equations(p):
    v, u = p
    return (v**4 - u**4 + 4*(c**2/3 + phi*(1+e))*v**2 - 4*(c**2/3+phi*(1-e))*u**2 + 16/3*e*phi*(phi-c**2), 1/(1+e)*v**3 - 1/(1-e)*u**3 + (3*phi+c**2/(1+e))*v-(3*phi+c**2/(1-e))*u)

v0 = np.sqrt(phi)*(1+e)
u0 = np.sqrt(phi)*(1-e)

print("Pericenter: ", a*(1-e))
print("Newtonian velocities: ", v0, "; ", u0)
vv, uu = fsolve(equations, (v0, u0))
print("Post-Newtonian velocities: ", vv, "; ", uu)