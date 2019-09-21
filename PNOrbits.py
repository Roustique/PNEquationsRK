#! /usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use("pgf")
#mpl.pyplot.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'], 'size': 16})
#mpl.pyplot.rc('font', family='Noto Serif')
pgf_with_custom_preamble = {
    "font.family": "sans-serif",
    "font.serif": [],
    "font.sans-serif": [],
    "text.usetex": True,
    "pgf.rcfonts": True,
    "pgf.preamble": [
#         "\\usepackage{units}",
#         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",
         r"\setmathfont{Latin Modern Math}[Scale=1.15]",
#        r"\setmainfont{lmodern}", # serif font via preamble
         ]
}
mpl.rcParams.update(pgf_with_custom_preamble)
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy import odr
from scipy.signal import argrelextrema
import datetime
import os

mpl.style.use('classic')
plt.rc('font', family='CMU Bright')
#plt.rc('text.latex',preamble=r'\usepackage[T2A]{fontenc}')
#plt.rc('text.latex',preamble=r'\usepackage[utf8]{inputenc}')
#plt.rc('text.latex',preamble=r'\usepackage[russian, english]{babel}')

pi = np.arctan(1.0) * 4.0
G = 1.0
timescale = 21.08
lengthscale = 6.32e+6


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
        self.timeweb = np.arange(self.size()) * self.dt

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
        return om

    def domega_turns(self):
        return self.omega()[1:] - self.omega()[:(np.size(self.omega()) - 1)]

    def domega_mean(self):
        return np.mean(self.domega_turns())

    def precession_turns(self):
        domegasum = np.empty(len(self.domega_turns()))
        for i in np.arange(len(self.domega_turns())):
            domegasum[i] = np.sum(self.domega_turns()[:i+1])
        return (domegasum / (self.timeweb[self.periloc()[1:]])) / timescale

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


def createoutputdir():
    outputdir = './output/' + datetime.datetime.now().strftime("%Y.%m.%d-%H:%M:%S") + '/'
    os.mkdir(outputdir)
    os.mkdir(outputdir+'anim/')
    os.system('cp ./input/data.dat ' + outputdir + 'data.dat')
    return outputdir

tstr = 'Time, '
cstr = 'Coordinate '
vstr = 'Velocity '
xstr = '$x$, '
ystr = '$y$, '
vxstr = '$v_x$, '
vystr = '$v_y$, '
lunit = 'km'
tunit = 'sec'
vunit = 'c'
Enstr = r'$\frac{E}{E_0}-1$'
Lnstr = r'$\frac{L}{L_0}-1$'
Orbitstr = 'Trajectory of system '
Depstr = 'Function '
xtstr = '$x(t)$ '
ytstr = '$y(t)$ '
vxtstr = '$v_x(t)$ '
vytstr = '$v_y(t)$ '
Eerr = 'Energy error '
Lerr = 'Angular moment error '


def outputIntegraline(Traj: Integraline, timescale, lengthscale, casename, casenameabbr, outputdir):
    casenamepdf = casenameabbr + '.pdf'

    print(casenameabbr+' - Drawing orbit trajectory...')
    fig, ax = plt.subplots()
    ax.scatter(0, 0, color='k')
    ax.plot(Traj.x * lengthscale, Traj.y * lengthscale, color='k')
    ax.axis('equal')
    ax.grid()
    ax.set(xlabel=cstr + xstr + lunit, ylabel=cstr + ystr + lunit, title=Orbitstr + casename)
    ax.margins(0.05)
    fig.savefig(outputdir + 'Orbit' + casenamepdf)

    print(casenameabbr+' - Drawing x(t)...')
    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.x * lengthscale, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + xstr + lunit, title=Depstr + xtstr + casename)
    fig.savefig(outputdir + 'x-t' + casenamepdf)

    print(casenameabbr+' - Drawing y(t)...')
    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.y * lengthscale, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + ystr + lunit, title=Depstr + ytstr + casename)
    fig.savefig(outputdir + 'y-t' + casenamepdf)

    print(casenameabbr+' - Drawing v_x(t)...')
    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.vx, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vxstr + vunit, title=Depstr + vxtstr + casename)
    fig.savefig(outputdir + 'vx-t' + casenamepdf)

    print(casenameabbr+' - Drawing v_y(t)...')
    fig, ax = plt.subplots()
    ax.plot(Traj.timeweb * timescale, Traj.vy, color='k')
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vystr + vunit, title=Depstr + vytstr + casename)
    fig.savefig(outputdir + 'vy-t' + casenamepdf)

    print(casenameabbr+' - Drawing Energy and angular moment plots...')
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
    
    print(casenameabbr+' - Drawing orbit trajectories...')
    fig, ax = plt.subplots()
    ax.scatter(0, 0, color='k')
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].x * lengthscale, Trajs[i].y * lengthscale, color=col[i])
    ax.axis('equal')
    ax.grid()
    ax.set(xlabel=cstr + xstr + lunit, ylabel=cstr + ystr + lunit, title=Orbitstr + casename)
    ax.margins(0.05)
    fig.savefig(outputdir + 'Orbit' + casenamepdf)

    print(casenameabbr+' - Drawing x(t)...')
    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].x * lengthscale, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + xstr + lunit, title=Depstr + xtstr + casename)
    fig.savefig(outputdir + 'x-t' + casenamepdf)

    print(casenameabbr+' - Drawing y(t)...')
    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].y * lengthscale, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=cstr + ystr + lunit, title=Depstr + ytstr + casename)
    fig.savefig(outputdir + 'y-t' + casenamepdf)

    print(casenameabbr+' - Drawing v_x(t)...')
    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].vx, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vxstr + vunit, title=Depstr + vxtstr + casename)
    fig.savefig(outputdir + 'vx-t' + casenamepdf)

    print(casenameabbr+' - Drawing v_y(t)...')
    fig, ax = plt.subplots()
    for i in np.arange(len(Trajs)):
        ax.plot(Trajs[i].timeweb * timescale, Trajs[i].vy, color=col[i])
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel=vstr + vystr + vunit, title=Depstr + vytstr + casename)
    fig.savefig(outputdir + 'vy-t' + casenamepdf)

    print(casenameabbr+' - Drawing Energies and angular moments plots...')
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
    print(casenameabbr+' - Drawing pericenter argument plot...')
    fig, ax = plt.subplots()
    ax.scatter(Traj.periloc0() * Traj.dt * timescale, Traj.omega()-pi/2, color="k")
    ax.plot(Traj.periloc0() * Traj.dt * timescale, Traj.omega()-pi/2, color="k")
    ax.grid()
    ax.margins(0.05)
    ax.set(xlabel=tstr + tunit, ylabel='Pericenter argument, radians', title='Pericenter shift')
    fig.savefig(outputdir + 'domega' + casenamepdf)

    print(casenameabbr+' - Drawing major semi-axis and eccentricity plot...')
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.scatter((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.a_turns() * lengthscale, color="k")
    ax1.plot((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.a_turns() * lengthscale, color="k")
    ax1.grid()
    ax1.margins(0.05)
    ax1.set(xlabel='Orbit turns', ylabel='Major semi-axis, km', title='Major semi-axis error')
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2 = plt.subplot(212)
    ax2.scatter((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.ecc_turns(), color="k")
    ax2.plot((np.arange(Traj.number_of_full_turns()) + 1).astype(int), Traj.ecc_turns(), color="k")
    ax2.grid()
    ax2.margins(0.05)
    ax2.set(xlabel='Orbit turns', ylabel='Eccentricity', title='Eccentricity error')
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    fig.savefig(outputdir + 'Keplerpar' + casenamepdf)


def whitemesh(colin: tuple, segnumb, shade):
    colar = np.array(colin)
    newcol = 1.0 - (1.0 - colar)/segnumb*shade
    return tuple(newcol)
    
def colmesh(col1: tuple, col2: tuple, segnumb, shade):
    col1ar = np.array(col1)
    col2ar = np.array(col2)
    newcol = col2ar - (col2ar - col1ar)/segnumb*shade
    return tuple(newcol)
    
def linethickness(segnumb, k):
    return 2.0 - 1.0*k/segnumb
    
def ellipsedots(a, e, omega):
    phi = np.arange(1000)/1000*2*pi
    r = a * (1 - e ** 2)/(1 + e * np.cos(phi))
    x = r * np.cos(phi + omega)
    y = r * np.sin(phi + omega)
    return np.array([x, y])
    
def outputframes1(Trajs: list, timescale, lengthscale, durationsec, tracelength, segnumb, outputdir):
    numberofframes = int(durationsec*24)
    col = [(0.25, 0.25, 0.25), (0.5, 0.25, 0), (0, 0.25, 0.5), (0, 0.5, 0.25)]
    colbright = [(1.0, 1.0, 1.0), (1.0, 0.5, 0), (0, 0.5, 1.0), (0, 1.0, 0.5)]
    whitecol=np.array([1.0, 1.0, 1.0])
    maxy=0
    maxx=0
    miny=0
    minx=0
    for j in np.arange(len(Trajs)):
        maxy = max(np.max(Trajs[j].y),maxy)
        maxx = max(np.max(Trajs[j].x),maxx)
        miny = min(np.min(Trajs[j].y),miny)
        minx = min(np.min(Trajs[j].x),minx)
    imsize = np.array([maxx-minx, maxy-miny])
    if (imsize[1]*4/3 > imsize[0]):
        imsize[0] = imsize[1]*4.0/3.0
    else:
        imsize[1] = imsize[0]*3.0/4.0
    imsize = imsize * 1.08
    centerpoint = np.array([maxx + minx, maxy + miny])/2
    ylimits = [(centerpoint[1]-imsize[1]/2)*lengthscale, (centerpoint[1]+imsize[1]/2)*lengthscale]
    xlimits = [(centerpoint[0]-imsize[0]/2)*lengthscale, (centerpoint[0]+imsize[0]/2)*lengthscale]
    for i in np.arange(numberofframes):
        fig, ax = plt.subplots()
        ax.axis('off')
        ax.set_aspect(aspect=1)
        fig.set_facecolor('k')
        ax.set_ylim(ylimits)
        ax.set_xlim(xlimits)
        ax.margins(0.01)
        ax.scatter([0], [0], color=(0.5, 0.25, 0), s=90, linewidths=0, zorder = 1)        
        ax.scatter([0], [0], color=(1.0, 0.5, 0), s=75, linewidths=0, zorder = 2)        
        ax.scatter([0], [0], color=(1.0, 1.0, 1.0), s=55, linewidths=0, zorder = 3)        
        ax.scatter([0], [0], color=(1.0, 0.5, 0), s=40, linewidths=0, zorder = 4)        
        ax.scatter([0], [0], color=(0, 0, 0), s=25, linewidths=0, zorder = 5)
        for j in np.arange(len(Trajs)):
            ppf = int(np.floor(Trajs[j].size()/numberofframes))  #points per frame
            cp = ppf*i   #current point
            ax.plot(Trajs[j].x * lengthscale, Trajs[j].y * lengthscale, color=col[j])
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = (linethickness(segnumb, k)-1)*4.5, color = colmesh((0, 0, 0), col[j], segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s=60, color = col[j], linewidths=0, zorder=96)
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = (linethickness(segnumb, k)-1)*3, color = colmesh(col[j], colbright[j], segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s = 30, color = colbright[j], linewidths=0, zorder=97)
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = linethickness(segnumb, k), color = colmesh(col[j], (1.0, 1.0, 1.0), segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s = 18, color = (1.0, 1.0, 1.0), linewidths=0, zorder=98)
        fig.savefig(outputdir+'anim/'+str(i)+'.png', facecolor='k')
        plt.close(fig)
        
def outputframes2(Trajs: list, timescale, lengthscale, durationsec, tracelength, segnumb, outputdir):
    numberofframes = int(durationsec*24)
    col = [(0.25, 0.25, 0.25), (0.5, 0.25, 0), (0, 0.25, 0.5), (0, 0.5, 0.25)]
    colbright = [(1.0, 1.0, 1.0), (1.0, 0.5, 0), (0, 0.5, 1.0), (0, 1.0, 0.5)]
    whitecol=np.array([1.0, 1.0, 1.0])
    a = np.empty(len(Trajs))
    e = np.empty(len(Trajs))
    prec = np.empty(len(Trajs))
    dt = np.empty(len(Trajs))
    for j in np.arange(len(Trajs)):
        a[j] = Trajs[j].a_mean()
        e[j] = Trajs[j].ecc_mean()
        prec[j] = Trajs[j].precession_mean()
        dt[j] = Trajs[j].dt
    maxy=0
    maxx=0
    miny=0
    minx=0
    for j in np.arange(len(Trajs)):
        maxy = max(np.max(Trajs[j].y),maxy)
        maxx = max(np.max(Trajs[j].x),maxx)
        miny = min(np.min(Trajs[j].y),miny)
        minx = min(np.min(Trajs[j].x),minx)
    imsize = np.array([maxx-minx, maxy-miny])
    if (imsize[1]*4/3 > imsize[0]):
        imsize[0] = imsize[1]*4.0/3.0
    else:
        imsize[1] = imsize[0]*3.0/4.0
    imsize = imsize * 1.08
    centerpoint = np.array([maxx + minx, maxy + miny])/2
    ylimits = [(centerpoint[1]-imsize[1]/2)*lengthscale, (centerpoint[1]+imsize[1]/2)*lengthscale]
    xlimits = [(centerpoint[0]-imsize[0]/2)*lengthscale, (centerpoint[0]+imsize[0]/2)*lengthscale]
    for i in np.arange(numberofframes):
        fig, ax = plt.subplots()
        ax.axis('off')
        ax.set_aspect(aspect=1)
        fig.set_facecolor('k')
        ax.set_ylim(ylimits)
        ax.set_xlim(xlimits)
        ax.margins(0.01)
        ax.scatter([0], [0], color=(0.5, 0.25, 0), s=90, linewidths=0, zorder = 1)        
        ax.scatter([0], [0], color=(1.0, 0.5, 0), s=75, linewidths=0, zorder = 2)        
        ax.scatter([0], [0], color=(1.0, 1.0, 1.0), s=55, linewidths=0, zorder = 3)        
        ax.scatter([0], [0], color=(1.0, 0.5, 0), s=40, linewidths=0, zorder = 4)        
        ax.scatter([0], [0], color=(0, 0, 0), s=25, linewidths=0, zorder = 5)
        for j in np.arange(len(Trajs)):
            ppf = int(np.floor(Trajs[j].size()/numberofframes))  #points per frame
            cp = ppf*i   #current point
            orbit = ellipsedots(a[j], e[j], prec[j] * cp * timescale * dt[j] + pi/2)
            ax.plot(orbit[0] * lengthscale, orbit[1] * lengthscale, color=col[j])
            #ax.plot(Trajs[j].x * lengthscale, Trajs[j].y * lengthscale, color=col[j])
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = (linethickness(segnumb, k)-1)*4.5, color = colmesh((0, 0, 0), col[j], segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s=60, color = col[j], linewidths=0, zorder=54)
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = (linethickness(segnumb, k)-1)*3, color = colmesh(col[j], colbright[j], segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s = 30, color = colbright[j], linewidths=0, zorder=56)
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = linethickness(segnumb, k), color = colmesh(col[j], (1.0, 1.0, 1.0), segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s = 18, color = (1.0, 1.0, 1.0), linewidths=0, zorder=58)
        fig.savefig(outputdir+'anim/'+str(i)+'.png', facecolor='k')
        plt.close(fig)
        
def outputframes3(Trajs: list, timescale, lengthscale, durationsec, tracelength, segnumb, outputdir):
    numberofframes = int(durationsec*24)
    col = [(0.25, 0.25, 0.25), (0.0, 0.25, 0.5), (0.5, 0.25, 0), (0.25, 0, 0.5)]
    colbright = [(1.0, 1.0, 1.0), (0.0, 0.5, 1.0), (1.0, 0.5, 0), (0.5, 0, 1.0)]
    whitecol=np.array([1.0, 1.0, 1.0])
    maxy=0
    maxx=0
    miny=0
    minx=0
    for j in np.arange(len(Trajs)):
        maxy = max(np.max(Trajs[j].y),maxy)
        maxx = max(np.max(Trajs[j].x),maxx)
        miny = min(np.min(Trajs[j].y),miny)
        minx = min(np.min(Trajs[j].x),minx)
    imsize = np.array([maxx-minx, maxy-miny])
    if (imsize[1]*4/3 > imsize[0]):
        imsize[0] = imsize[1]*4.0/3.0
    else:
        imsize[1] = imsize[0]*3.0/4.0
    imsize = imsize * 1.08
    centerpoint = np.array([maxx + minx, maxy + miny])/2
    ylimits = [(centerpoint[1]-imsize[1]/2)*lengthscale, (centerpoint[1]+imsize[1]/2)*lengthscale]
    xlimits = [(centerpoint[0]-imsize[0]/2)*lengthscale, (centerpoint[0]+imsize[0]/2)*lengthscale]
    for i in np.arange(numberofframes):
        fig, ax = plt.subplots()
        ax.axis('off')
        ax.set_aspect(aspect=1)
        fig.set_facecolor('k')
        ax.set_ylim(ylimits)
        ax.set_xlim(xlimits)
        ax.margins(0.01)
        ax.scatter([0], [0], color=(0.5, 0.25, 0), s=90, linewidths=0, zorder = 1)        
        ax.scatter([0], [0], color=(1.0, 0.5, 0), s=75, linewidths=0, zorder = 2)        
        ax.scatter([0], [0], color=(1.0, 1.0, 1.0), s=55, linewidths=0, zorder = 3)        
        ax.scatter([0], [0], color=(1.0, 0.5, 0), s=40, linewidths=0, zorder = 4)        
        ax.scatter([0], [0], color=(0, 0, 0), s=25, linewidths=0, zorder = 5)
        for j in np.arange(len(Trajs)):
            ppf = int(np.floor(Trajs[j].size()/numberofframes))  #points per frame
            cp = ppf*i   #current point
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = (linethickness(segnumb, k)-1)*4, color = colmesh((0, 0, 0), col[j], segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s=60, color = col[j], linewidths=0, zorder=54)
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = (linethickness(segnumb, k)-1)*2, color = colmesh(col[j], colbright[j], segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s = 30, color = colbright[j], linewidths=0, zorder=56)
            for k in np.arange(segnumb):
                br = cp - ppf*(k+1)*tracelength  #bottom range
                br = int(br * (br > 0))
                ur = cp - ppf*k*tracelength      #upper range
                ur = int(ur * (ur > 0))
                ax.plot((Trajs[j].x)[br:ur] * lengthscale, (Trajs[j].y)[br:ur] * lengthscale, linewidth = linethickness(segnumb, k)*0.6, color = colmesh(col[j], (1.0, 1.0, 1.0), segnumb, k))
            ax.scatter(Trajs[j].x[cp] * lengthscale, Trajs[j].y[cp] * lengthscale, s = 18, color = (1.0, 1.0, 1.0), linewidths=0, zorder=58)
        fig.savefig(outputdir+'anim/'+str(i)+'.png', facecolor='k')
        plt.close(fig)
