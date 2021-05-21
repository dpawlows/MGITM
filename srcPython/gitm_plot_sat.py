#!/usr/bin/env python

from glob import glob
from datetime import datetime
from datetime import timedelta
from struct import unpack
import numpy as np
import matplotlib.pyplot as pp
import matplotlib.dates as dates
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
from gitm_routines import *
import re
import sys


files = glob("1D*")
filetypes = list(set([f[0:5] for f in files]))
filetypes.sort()

for i in range(len(filetypes)):
    print("{}: {}".format(i,filetypes[i]))

itype = int(input("Which filetype to plot:"))

filestoplot = [thisfile for thisfile in files if thisfile[0:5] == filetypes[itype]]
filestoplot.sort()
nfiles = len(filestoplot)

header = read_gitm_header(filestoplot)

for i in range(header['nVars']):
    print("{}: {}".format(i,header['vars'][i]))

pvar = int(input("Which variable to plot: "))
vars = [0,1,2,pvar]
time = []

for ifile in range(len(filestoplot)):
    data = read_gitm_one_file(filestoplot[ifile],vars)
    time.append(data["time"])
    if ifile == 0:
        altitude = data[2][0][0]/1000. #Convert to km
        plotdata = np.zeros((nfiles,data['nAlts']))

    plotdata[ifile,:] = data[15]

X,Y = np.meshgrid(time,altitude)
fig, ax = plt.subplots()
cont = ax.contourf(X,Y,plotdata.transpose(),cmap="jet",levels=30)
pp.ylabel('Altitude (km)')
# pp.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y %H:%M:%S'))
fig.autofmt_xdate()

cbar = pp.colorbar(cont)
cbar.set_label('{}'.format(header['vars'][15]))

pp.savefig('plot.png')
