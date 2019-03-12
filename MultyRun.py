#! /usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import scipy
from scipy import odr
from scipy.signal import argrelextrema
import math
import datetime
import os
import subprocess

#plt.rc('font', family='Computer Modern')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 16})
plt.rc('text', usetex=True)
plt.rc('text.latex',unicode=True)
plt.rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
plt.rc('text.latex',preamble=r'\usepackage[russian]{babel}')


inputfile = open('./input/dataMulty.dat', 'r')
inputfile.readline()
for line in inputfile:
    l = line.split(sep=", ")
    l[6]=l[6].replace("\n", "")
    
    outputfile = open('./input/data.dat', 'w')
    outputfile.write(l[0]+"\n")
    outputfile.write(l[1]+"\n")
    outputfile.write(l[2]+"\n")
    outputfile.write(l[3]+"\n")
    outputfile.write(l[4]+"\n")
    outputfile.write(l[5]+"\n")
    outputfile.write(l[6])
    outputfile.close()
    
    subprocess.call(["./main.out", "-N"])
    subprocess.call(["./main.out", "-P"])
    subprocess.call(["python3", "PlotDrawer.py"])