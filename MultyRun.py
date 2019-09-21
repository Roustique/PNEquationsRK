#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess

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
    
    print('Running main.out -NT ...')
    subprocess.call(["./main.out", "-NT"])
    print('done!'+'\n'+'Running main.out -GS ...')
    subprocess.call(["./main.out", "-GS"])
    print('done!'+'\n'+'Running main.out -GH ...')
    subprocess.call(["./main.out", "-GH"])
    print('done!'+'\n'+'Running main.out -F0 ...')
    subprocess.call(["./main.out", "-F0"])
    print('done!'+'\n'+'Running PlotDrawer.py ...')
    subprocess.call(["python3", "PlotDrawer.py"])
    print('done!')