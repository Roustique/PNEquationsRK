#! /usr/bin/env python
# -*- coding: utf-8 -*-

import PNOrbits

timescale = 21.08
lengthscale = 6.32e+6

inputfile=open('./input/data.dat')
dt = float(inputfile.readline())
inputfile.readline()
inputfile.readline()
inputfile.readline()
inputfile.readline()
inputfile.readline()
M = float(inputfile.readline())

TrajNT = PNOrbits.Integraline(dt, M)
TrajGS = PNOrbits.Integraline(dt, M)
TrajGH = PNOrbits.Integraline(dt, M)
TrajF0 = PNOrbits.Integraline(dt, M)

print('NT - Reading output result')
TrajNT.readtxt('./output/result-NT.dat')
print('GS - Reading output result')
TrajGS.readtxt('./output/result-GS.dat')
print('GH - Reading output result')
TrajGH.readtxt('./output/result-GH.dat')
print('F0 - Reading output result')
TrajF0.readtxt('./output/result-F0.dat')

outputdir = PNOrbits.createoutputdir()

PNOrbits.outputIntegraline(TrajNT, timescale, lengthscale, '(Newtonian)', 'NT', outputdir)
PNOrbits.outputIntegraline(TrajGS, timescale, lengthscale, '(GR, Schwarzschild coord.)', 'GS', outputdir)
PNOrbits.outputIntegraline(TrajGH, timescale, lengthscale, '(GR, harm. coord; FGT, $K=1/2$)', 'GH', outputdir)
PNOrbits.outputIntegraline(TrajF0, timescale, lengthscale, '(FGT, $K=0$)', 'F0', outputdir)
PNOrbits.outputIntegralines([TrajNT, TrajGS, TrajGH, TrajF0], timescale, lengthscale, '(All)', 'All', outputdir)
PNOrbits.outputPars(TrajGS, lengthscale, '(GR, Schwarzschild coord.)', 'GS', outputdir)
PNOrbits.outputPars(TrajGH, lengthscale, '(GR, harm. coord; FGT, $K=1/2$)', 'GH', outputdir)
PNOrbits.outputPars(TrajF0, lengthscale, '(FGT, $K=0$)', 'F0', outputdir)
