#! /usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.font_manager as mfm
import numpy as np
from scipy import odr
from scipy.signal import argrelextrema
import datetime
import os

plt.rc('font', family='Geometria')

pi = np.arctan(1.0) * 4.0
G = 1.47662e-1
timescale = 3.335640952e-5
lengthscale = 10.0


def ellipsepolar(b, phi):
    return b[0] * (1 - b[1] ** 2) / (1 + b[1] * np.cos(phi))


class Integraline:

    def __init__(self, dt, M):
        self.x = np.empty(0)
        self.y = np.empty(0)
        self.vx = np.empty(0)
        self.vy = np.empty(0)
        self.E = np.empty(0)
        self.L = np.empty(0)
        self.r = np.empty(0)
        self.phi = np.empty(0)
        self.vabs = np.empty(0)
        self.Enormed = np.empty(0)
        self.Lnormed = np.empty(0)
        self.timeweb = np.empty(0)
        self.dt = dt
        self.M = M


    def readinputpars(self, file=open('./input/data.dat')):
        self.dt = float(file.readline())
        file.readline()
        file.readline()
        file.readline()
        file.readline()
        file.readline()
        self.M = float(file.readline())

    def readtxt(self, filename):
        self.x = np.loadtxt(open(filename, 'r'), delimiter=", ", skiprows=1, usecols=(0,))
        self.y = np.loadtxt(open(filename, 'r'), delimiter=", ", skiprows=1, usecols=(1,))
        self.vx = np.loadtxt(open(filename, 'r'), delimiter=", ", skiprows=1, usecols=(2,))
        self.vy = np.loadtxt(open(filename, 'r'), delimiter=", ", skiprows=1, usecols=(3,))
        self.E = np.loadtxt(open(filename, 'r'), delimiter=", ", skiprows=1, usecols=(4,))
        self.L = np.loadtxt(open(filename, 'r'), delimiter=", ", skiprows=1, usecols=(5,))
        self.Enormed = self.E / self.E[0] - 1
        self.Lnormed = self.L / self.L[0] - 1
        self.r = np.sqrt(self.x ** 2 + self.y ** 2)
        self.phi = np.arctan2(self.y, self.x)
        for i in np.arange(1, self.size()):
            if 0 > self.phi[i] - self.phi[i - 1]:
                self.phi[i:] = self.phi[i:] + 2 * pi
        self.vabs = np.sqrt(self.vx ** 2 + self.vy ** 2)
        self.timeweb = np.arange(self.size()) * dt

    def size(self):
        return np.size(self.x)

    def periloc0(self):
        return np.insert(argrelextrema(self.r, np.less), 0, 0)

    def number_of_full_turns(self):
        return np.size(self.periloc0()) - 1

    def turn_length(self):
        n = np.size(self.periloc0())
        return self.periloc0()[1:] - self.periloc0()[:n-1]

    def turn_length_mean(self):
        return np.mean(self.turn_length()[np.arange(self.number_of_full_turns())])

    def periloc(self):
        if abs(self.turn_length_mean() - self.periloc0()[1]) < 2:
            return self.periloc0()
        else:
            return self.periloc0()[1:]

    def omega(self):
        om = np.arctan2(self.y[self.periloc()], self.x[self.periloc()])
        for i in np.arange(1, self.number_of_full_turns() + 1):
            if om[i] - om[i - 1] < 0:
                om[i:] = om[i:] + 2 * pi
        return om

    def domega_turns(self):
        return self.omega()[1:] - self.omega()[:(np.size(self.omega()) - 1)]

    def domega_mean(self):
        return np.mean(self.domega_turns())

    def precession_turns(self):
        return (self.domega_turns()[1:] / self.timeweb[self.periloc()[1:]]) / timescale

    def precession_mean(self):
        return np.mean(self.precession_turns())

    def phi_rotated(self):
        phir = np.empty(self.size())
        for i in np.arange(self.number_of_full_turns()):
            phir[self.periloc()[i]:self.periloc()[i + 1]] = (self.phi - (
                    self.domega_turns()[i] / (self.timeweb[self.periloc()[i + 1]] - self.timeweb[self.periloc()[i]])) * self.timeweb - self.omega()[
                                                                 0])[self.periloc()[i]:self.periloc()[i + 1]]
        phir[self.periloc()[self.number_of_full_turns()]:] = (self.phi - (
                self.domega_turns()[self.number_of_full_turns()-1] / (self.timeweb[
            self.periloc()[self.number_of_full_turns()]]-self.timeweb[self.periloc()[self.number_of_full_turns()-1]])) * self.timeweb - self.omega()[0])[
                                                             self.periloc()[self.number_of_full_turns()]:]
        return phir

    def fit_to_model(self):
        ellipsemodel = odr.Model(ellipsepolar)
        a = np.empty(self.number_of_full_turns())
        ecc = np.empty(self.number_of_full_turns())
        for i in np.arange(self.number_of_full_turns()):
            data_to_fit = odr.Data(self.phi_rotated()[self.periloc()[i]:self.periloc()[i + 1]],
                                   self.r[self.periloc()[i]:self.periloc()[i + 1]])
            peri0 = self.r[self.periloc()[0]]
            apo0 = self.r[int(self.periloc()[0] + self.turn_length()[0] / 2)]
            a0 = (peri0 + apo0) / 2
            ecc0 = (1 - peri0 / apo0) / (1 + peri0 / apo0)
            job = odr.ODR(data_to_fit, ellipsemodel, beta0=np.array([a0, ecc0]))
            results = job.run()
            a[i], ecc[i] = results.beta
        return a, ecc

    def a_turns(self):
        return self.fit_to_model()[0]

    def ecc_turns(self):
        return self.fit_to_model()[1]

    def a_mean(self):
        return np.mean(self.a_turns())

    def ecc_mean(self):
        return np.mean(self.ecc_turns())

    def a_sd(self):
        return np.sqrt(np.sum((self.a_turns() - self.a_mean()) ** 2))

    def ecc_sd(self):
        return np.sqrt(np.sum((self.ecc_turns() - self.ecc_mean()) ** 2))


inputfile=open('./input/data.dat')
dt = float(inputfile.readline())
inputfile.readline()
inputfile.readline()
inputfile.readline()
inputfile.readline()
inputfile.readline()
M = float(inputfile.readline())

TrajN = Integraline(dt, M)
TrajP = Integraline(dt, M)
TrajN.readtxt('./output/resultN.dat')
TrajP.readtxt('./output/resultP.dat')

outputdir = './output/' + datetime.datetime.now().strftime("%Y.%m.%d-%H:%M:%S") + '/'
os.mkdir(outputdir)
os.system('cp ./input/data.dat ' + outputdir + 'data.dat')

tstr = 'Время, '
cstr = 'Координата '
vstr = 'Скорость '
xstr = '$x$, '
ystr = '$y$, '
vxstr = '$v_x$, '
vystr = 'v, '
lunit = 'км'
tunit = 'сек'
vunit = 'c'
Enstr = r'$\frac{E}{E_0}-1$'
Lnstr = r'$\frac{L}{L_0}-1$'
Orbitstr = 'Траектория системы '
Depstr = 'Зависимость '
xtstr = '$x(t)$ '
ytstr = '$y(t)$ '
vxtstr = '$v_x(t)$ '
vytstr = '$v_y(t)$ '
Eerr = 'Ошибки энергий '
Lerr = 'Ошибки угловых моментов '


def outputIntegraline(Traj: Integraline, timescale, lengthscale, casename, casenameabbr, outputdir):
    casenamepdf = casenameabbr + '.pdf'

    fig, ax = plt.subplots()
    ax.scatter(0, 0, color='k')
    ax.plot(Traj.x * lengthscale, Traj.y * lengthscale, color='k')
    ax.axis('equal')
    ax.grid()
    ax.set(xlabel=cstr + xstr + lunit, ylabel=cstr + ystr + lunit, title=Orbitstr + casename)
    ax.margins(0.05)
    fig.savefig(outputdir + 'Orbit' + casenamepdf)

    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.x * lengthscale, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + xstr + lunit, title=Depstr + xtstr + casename)
    fig.savefig(outputdir + 'x-t' + casenamepdf)

    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.y * lengthscale, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + ystr + lunit, title=Depstr + ytstr + casename)
    fig.savefig(outputdir + 'y-t' + casenamepdf)

    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.vx, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vxstr + vunit, title=Depstr + vxtstr + casename)
    fig.savefig(outputdir + 'vx-t' + casenamepdf)

    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.vy, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vystr + vunit, title=Depstr + vytstr + casename)
    fig.savefig(outputdir + 'vy-t' + casenamepdf)

    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.plot(Traj.timeweb * timescale, Traj.Enormed, color='k')
    ax1.grid()
    ax1.margins(0.05)
    ax1.set(xlabel=tstr + tunit, ylabel=Enstr, title=Eerr + casename)
    ax2 = plt.subplot(212)
    ax2.plot(Traj.timeweb * timescale, Traj.Lnormed, color='k')
    ax2.grid()
    ax2.margins(0.05)
    ax2.set(xlabel=tstr + tunit, ylabel=Lnstr, title=Lerr + casename)
    plt.tight_layout()
    fig.savefig(outputdir + 'Errors' + casenamepdf)


def outputIntegralines(Trajs: list, timescale, lengthscale, casename, casenameabbr, outputdir):
    """

    :type Trajs: list[Integraline]
    """
    casenamepdf = casenameabbr + '.pdf'
    col = ['k', 'r', 'b', 'g', 'o']
    col = col[:len(Trajs)]
    fig, ax = plt.subplots()
    ax.scatter(0, 0, color='k')
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].x * lengthscale, Trajs[i].y * lengthscale, color=col[i])
    ax.axis('equal')
    ax.grid()
    ax.set(xlabel=cstr + xstr + lunit, ylabel=cstr + ystr + lunit, title=Orbitstr + casename)
    ax.margins(0.05)
    fig.savefig(outputdir + 'Orbit' + casenamepdf)

    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].x * lengthscale, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + xstr + lunit, title=Depstr + xtstr + casename)
    fig.savefig(outputdir + 'x-t' + casenamepdf)

    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].y * lengthscale, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + ystr + lunit, title=Depstr + ytstr + casename)
    fig.savefig(outputdir + 'y-t' + casenamepdf)

    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].vx, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vxstr + vunit, title=Depstr + vxtstr + casename)
    fig.savefig(outputdir + 'vx-t' + casenamepdf)

    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].vy, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vystr + vunit, title=Depstr + vytstr + casename)
    fig.savefig(outputdir + 'vy-t' + casenamepdf)

    fig = plt.figure()
    ax1 = plt.subplot(211)
    for i in np.arange(len(Trajs)):
        ax1.plot(Trajs[i].timeweb * timescale, Trajs[i].Enormed, color=col[i])
    ax1.grid()
    ax1.margins(0.05)
    ax1.set(xlabel=tstr + tunit, ylabel=Enstr, title=Eerr + casename)
    ax2 = plt.subplot(212)
    for i in np.arange(len(Trajs)):
        ax2.plot(Trajs[i].timeweb * timescale, Trajs[i].Lnormed, color=col[i])
    ax2.grid()
    ax2.margins(0.05)
    ax2.set(xlabel=tstr + tunit, ylabel=Lnstr, title=Lerr + casename)
    plt.tight_layout()
    fig.savefig(outputdir + 'Errors' + casenamepdf)


def outputPars(Traj: Integraline, lengthscale, casename, casenameabbr, outputdir):
    casenamepdf = casenameabbr + '.pdf'
    fig, ax = plt.subplots()
    ax.scatter(Traj.periloc0() * Traj.dt * timescale, Traj.omega()-pi/2, color="k")
    ax.plot(Traj.periloc0() * Traj.dt * timescale, Traj.omega()-pi/2, color="k")
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel='Аргумент перицентра, радианы', title='Изменение аргумента перицентра')
    fig.savefig(outputdir + 'domega' + casenamepdf)

    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.scatter((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.a_turns() * lengthscale, color="k")
    ax1.plot((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.a_turns() * lengthscale, color="k")
    ax1.grid()
    ax1.margins(0.05)
    ax1.set(xlabel='Обороты по орбите', ylabel='Большая полуось, км', title='Изменение большой полуоси')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2 = plt.subplot(212)
    ax2.scatter((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.ecc_turns(), color="k")
    ax2.plot((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.ecc_turns(), color="k")
    ax2.grid()
    ax2.margins(0.05)
    ax2.set(xlabel='Обороты по орбите', ylabel='Эксцентриситет', title='Изменение эксцентриситета')
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    fig.savefig(outputdir + 'Keplerpar' + casenamepdf)


outputIntegraline(TrajN, timescale, lengthscale, '(Н. случай)', 'N', outputdir)
outputIntegraline(TrajP, timescale, lengthscale, '(П-Н. случай)', 'P', outputdir)
outputIntegralines([TrajN, TrajP], timescale, lengthscale, '(Н. и П-Н. случаи)', 'NP', outputdir)
outputPars(TrajP, lengthscale, '(П-Н. случай)', 'P', outputdir)
