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

#plt.rc('font', family='Computer Modern')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 16})
plt.rc('text', usetex=True)
plt.rc('text.latex',unicode=True)
plt.rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
plt.rc('text.latex',preamble=r'\usepackage[russian]{babel}')


def ellipsepolar(B, phi):
    return B[0]*(1-B[1]**2)/(1+B[1]*np.cos(phi))


xN = np.loadtxt(open('./output/resultN.dat', 'r'), delimiter=', ', skiprows=1, usecols=(0,))
yN = np.loadtxt(open('./output/resultN.dat', 'r'), delimiter=', ', skiprows=1, usecols=(1,))
vxN = np.loadtxt(open('./output/resultN.dat', 'r'), delimiter=', ', skiprows=1, usecols=(2,))
vyN = np.loadtxt(open('./output/resultN.dat', 'r'), delimiter=', ', skiprows=1, usecols=(3,))
EN = np.loadtxt(open('./output/resultN.dat', 'r'), delimiter=', ', skiprows=1, usecols=(4,))
LN = np.loadtxt(open('./output/resultN.dat', 'r'), delimiter=', ', skiprows=1, usecols=(5,))
xP = np.loadtxt(open('./output/resultP.dat', 'r'), delimiter=', ', skiprows=1, usecols=(0,))
yP = np.loadtxt(open('./output/resultP.dat', 'r'), delimiter=', ', skiprows=1, usecols=(1,))
vxP = np.loadtxt(open('./output/resultP.dat', 'r'), delimiter=', ', skiprows=1, usecols=(2,))
vyP = np.loadtxt(open('./output/resultP.dat', 'r'), delimiter=', ', skiprows=1, usecols=(3,))
EP = np.loadtxt(open('./output/resultP.dat', 'r'), delimiter=', ', skiprows=1, usecols=(4,))
LP = np.loadtxt(open('./output/resultP.dat', 'r'), delimiter=', ', skiprows=1, usecols=(5,))

inputfile = open('./input/data.dat')
dt = float(inputfile.readline())
n = int(inputfile.readline())
inputfile.readline()
inputfile.readline()
inputfile.readline()
inputfile.readline()
M = float(inputfile.readline())
G = 1.47662e-1
pi = np.arctan(1.0)*4.0
timeweb = np.arange(n+1)*dt
rN = np.sqrt(xN**2+yN**2)
vabsN = np.sqrt(vxN**2+vyN**2)
aN = rN[0]*G*M/(2*G*M-vabsN[0]**2*rN[0])    #Semi-major axis
omegaN = math.atan2(yN[np.argmin(rN)], xN[np.argmin(rN)])  #Pericenter argument
eccN = np.max(np.roots(np.array([1, (xN[0]*np.cos(omegaN)+yN[0]*np.sin(omegaN))/aN, rN[0]/aN-1])))  #eccentricity

rP = np.sqrt(xP**2+yP**2)
periarg = np.insert(argrelextrema(rP, np.less), 0, 0)
numberofturns = np.size(periarg)-1
domega = (np.arctan2(yP[periarg], xP[periarg])-omegaN)
for i in np.arange(1,numberofturns+1):
    if (domega[i]-domega[i-1]<0):
        domega[i:]=domega[i:]+2*pi
phiP = np.arctan2(yP, xP)
phiProtated = np.empty(n+1)
for i in np.arange(numberofturns):
    phiProtated[periarg[i]:periarg[i+1]] = (phiP-(domega[i+1]/timeweb[periarg[i+1]])*timeweb-omegaN)[periarg[i]:periarg[i+1]]
phiProtated[periarg[numberofturns]:] = (phiP-(domega[numberofturns]/timeweb[periarg[numberofturns]])*timeweb-omegaN)[periarg[numberofturns]:]
aP = np.empty(numberofturns)
eccP = np.empty(numberofturns)

for i in np.arange(numberofturns):
    ellipsemodel = odr.Model(ellipsepolar)
    data_to_fit = odr.Data(phiProtated[periarg[i]:periarg[i+1]], rP[periarg[i]:periarg[i+1]])
    job = odr.ODR(data_to_fit, ellipsemodel, beta0=np.array([aN, eccN]))
    results = job.run()
    aP[i], eccP[i] = results.beta
    

#rNewtAnalyt = a*(1-ecc**2)/(1+ecc*(xN*np.cos(omega)+yN*np.sin(omega))/rN)
#xNewtAnalyt = xN/rN*rNewtAnalyt
#yNewtAnalyt = yN/rN*rNewtAnalyt

timescale = 3.335640952e-5
lengthscale = 10.0

EnormedN = EN/EN[0]-1
LnormedN = LN/LN[0]-1
EnormedP = EP/EP[0]-1
LnormedP = LP/LP[0]-1

outputdir = './output/'+datetime.datetime.now().strftime("%Y.%m.%d-%H:%M:%S")+'/'
os.mkdir(outputdir)
os.system('cp ./input/data.dat '+outputdir+'data.dat')

fig, ax = plt.subplots()
ax.scatter(0, 0, color='k')
ax.plot(xN*lengthscale, yN*lengthscale, color='k')
ax.axis('equal')
ax.grid()
ax.set(xlabel='Координата $x$, км', ylabel='Координата $y$, км', title='Орбита системы (Н. случай)')
ax.margins(0.05)
fig.savefig(outputdir+'OrbitN.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, xN*lengthscale, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Координата $x$, км', title='зависимость $x(t)$ (Н. случай)')
fig.savefig(outputdir+'x-tN.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, yN*lengthscale, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Координата $y$, км', title='зависимость $y(t)$ (Н. случай)')
fig.savefig(outputdir+'y-tN.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, vxN, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Скорость $v_x$, $c$', title='зависимость $v_x(t)$ (Н. случай)')
fig.savefig(outputdir+'vx-tN.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, vyN, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Скорость $v_y$, $c$', title='зависимость $v_y(t)$ (Н. случай)')
fig.savefig(outputdir+'vy-tN.png')

fig = plt.figure(6)
ax1 = plt.subplot(211)
ax1.plot(timeweb*timescale, EnormedN, color='k')
ax1.grid()
ax1.margins(0.05)
ax1.set(xlabel='Время, сек', ylabel=r'$\frac{E}{E_0}-1$', title='Ошибки энергии (Н. случай)')
#fig.savefig(outputdir+'EnergyError.png')
#fig, ax = plt.subplots()
ax2 = plt.subplot(212)
ax2.plot(timeweb*timescale, LnormedN, color='k')
ax2.grid()
ax2.margins(0.05)
ax2.set(xlabel='Время, сек', ylabel=r'$\frac{L}{L_0}-1$', title='Ошибки углового момента (Н. случай)')
plt.tight_layout()
fig.savefig(outputdir+'ErrorsN.png')


fig, ax = plt.subplots()
ax.scatter(0, 0, color='k')
ax.plot(xP*lengthscale, yP*lengthscale, color='k')
ax.axis('equal')
ax.grid()
ax.set(xlabel='Координата $x$, км', ylabel='Координата $y$, км', title='Орбита системы (П-н. случай)')
ax.margins(0.05)
fig.savefig(outputdir+'OrbitP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, xP*lengthscale, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Координата $x$, км', title='зависимость $x(t)$ (П-н. случай)')
fig.savefig(outputdir+'x-tP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, yP*lengthscale, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Координата $y$, км', title='зависимость $y(t)$ (П-н. случай)')
fig.savefig(outputdir+'y-tP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, vxP, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Скорость $v_x$, $c$', title='зависимость $v_x(t)$ (П-н. случай)')
fig.savefig(outputdir+'vx-tP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, vyP, color='k')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Скорость $v_y$, $c$', title='зависимость $v_y(t)$ (П-н. случай)')
fig.savefig(outputdir+'vy-tP.png')

fig = plt.figure(12)
ax1 = plt.subplot(211)
ax1.plot(timeweb*timescale, EnormedP, color='k')
ax1.grid()
ax1.margins(0.05)
ax1.set(xlabel='Время, сек', ylabel=r'$\frac{E}{E_0}-1$', title='Ошибки энергии (П-н. случай)')
#fig.savefig(outputdir+'EnergyError.png')
#fig, ax = plt.subplots()
ax2 = plt.subplot(212)
ax2.plot(timeweb*timescale, LnormedP, color='k')
ax2.grid()
ax2.margins(0.05)
ax2.set(xlabel='Время, сек', ylabel=r'$\frac{L}{L_0}-1$', title='Ошибки углового момента (П-н. случай)')
plt.tight_layout()
fig.savefig(outputdir+'ErrorsP.png')


fig, ax = plt.subplots()
ax.scatter(0, 0, color='k')
ax.plot(xN*lengthscale, yN*lengthscale, color='k')
ax.plot(xP*lengthscale, yP*lengthscale, color='r')
ax.axis('equal')
ax.grid()
ax.set(xlabel='Координата $x$, км', ylabel='Координата $y$, км', title='Орбиты системы')
ax.margins(0.05)
fig.savefig(outputdir+'OrbitNP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, xN*lengthscale, color='k')
ax.plot(timeweb*timescale, xP*lengthscale, color='r')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Координата $x$, км', title='зависимости $x(t)$')
fig.savefig(outputdir+'x-tNP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, yN*lengthscale, color='k')
ax.plot(timeweb*timescale, yP*lengthscale, color='r')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Координата $y$, км', title='зависимости $y(t)$')
fig.savefig(outputdir+'y-tNP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, vxN, color='k')
ax.plot(timeweb*timescale, vxP, color='r')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Скорость $v_x$, $c$', title='зависимости $v_x(t)$')
fig.savefig(outputdir+'vx-tNP.png')

fig, ax = plt.subplots()
ax.plot(timeweb*timescale, vyN, color='k')
ax.plot(timeweb*timescale, vyP, color='r')
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Скорость $v_y$, $c$', title='зависимости $v_y(t)$')
fig.savefig(outputdir+'vy-tNP.png')

fig = plt.figure(18)
ax1 = plt.subplot(211)
ax1.plot(timeweb*timescale, EnormedN, color='k')
ax1.plot(timeweb*timescale, EnormedP, color='r')
ax1.grid()
ax1.margins(0.05)
ax1.set(xlabel='Время, сек', ylabel=r'$\frac{E}{E_0}-1$', title='Ошибки энергий')
#fig.savefig(outputdir+'EnergyError.png')
#fig, ax = plt.subplots()
ax2 = plt.subplot(212)
ax2.plot(timeweb*timescale, LnormedN, color='k')
ax2.plot(timeweb*timescale, LnormedP, color='r')
ax2.grid()
ax2.margins(0.05)
ax2.set(xlabel='Время, сек', ylabel=r'$\frac{L}{L_0}-1$', title='Ошибки угловых моментов')
plt.tight_layout()
fig.savefig(outputdir+'ErrorsNP.png')


fig, ax = plt.subplots()
ax.scatter(periarg*dt*timescale, domega, color="k")
ax.plot(periarg*dt*timescale, domega, color="k")
ax.grid()
ax.margins(0.05)
ax.set(xlabel='Время, сек', ylabel='Аргумент перицентра, радианы', title='Изменение аргумента перицентра')
fig.savefig(outputdir+'domega.png')

fig = plt.figure(20)
ax1 = plt.subplot(211)
ax1.scatter((np.arange(numberofturns)+1).astype(int), aP*lengthscale, color="k")
ax1.plot((np.arange(numberofturns)+1).astype(int), aP*lengthscale, color="k")
ax1.grid()
ax1.margins(0.05)
ax1.set(xlabel='Обороты по орбите', ylabel='Большая полуось, км', title='Изменение большой полуоси')
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2 = plt.subplot(212)
ax2.scatter((np.arange(numberofturns)+1).astype(int), eccP, color="k")
ax2.plot((np.arange(numberofturns)+1).astype(int), eccP, color="k")
ax2.grid()
ax2.margins(0.05)
ax2.set(xlabel='Обороты по орбите', ylabel='Эксцентриситет', title='Изменение эксцентриситета')
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
fig.savefig(outputdir+'Keplerparchange.png')