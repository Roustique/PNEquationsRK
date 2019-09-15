#! /usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.font_manager as mfm
import matplotlib.style
import numpy as np
from scipy import odr
from scipy.signal import argrelextrema
import datetime
import subprocess
import os
import PNOrbits

matplotlib.style.use('classic')
plt.rc('font', family='Geometria')

timescale = 3.335640952e-5
lengthscale = 10.0

inputfile=open('./input/data.dat')
dt = float(inputfile.readline())
inputfile.readline()
inputfile.readline()
inputfile.readline()
inputfile.readline()
inputfile.readline()
M = float(inputfile.readline())

TrajN = PNOrbits.Integraline(dt, M)
TrajP = PNOrbits.Integraline(dt, M)
print('N - Reading output result')
TrajN.readtxt('./output/resultN.dat')
print('P - Reading output result')
TrajP.readtxt('./output/resultP.dat')

outputdir = PNOrbits.createoutputdir()

PNOrbits.outputframes3([TrajN, TrajP], timescale, lengthscale, 15, 0.3, 10, outputdir)
#subprocess.call(['bash', '-c', 'convert -delay 4.1667 -loop 0 '+outputdir+'anim/{0..'+str(6*24-1)+'}.png '+outputdir+'Output.gif'])
