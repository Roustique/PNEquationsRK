#! /usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.style
import matplotlib.font_manager as mfm
import numpy as np
from scipy import odr
from scipy.signal import argrelextrema
import datetime
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

PNOrbits.outputIntegraline(TrajN, timescale, lengthscale, '(Н. случай)', 'N', outputdir)
PNOrbits.outputIntegraline(TrajP, timescale, lengthscale, '(П-Н. случай)', 'P', outputdir)
PNOrbits.outputIntegralines([TrajN, TrajP], timescale, lengthscale, '(Н. и П-Н. случаи)', 'NP', outputdir)
PNOrbits.outputPars(TrajP, lengthscale, '(П-Н. случай)', 'P', outputdir)
