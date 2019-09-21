#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
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

PNOrbits.outputframes3([TrajNT, TrajGS, TrajGH, TrajF0], timescale, lengthscale, 12, 2.5, 10, outputdir)
subprocess.call(['bash', '-c', 'ffmpeg -framerate 24 -i '+outputdir+'anim/%00d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p '+outputdir+'video_output.mp4'])
