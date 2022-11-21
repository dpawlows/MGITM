#!/usr/bin/env python
#Plots a altitude and time for a single location
from glob import glob
from datetime import datetime
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as pp
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from gitm_routines import *
import re
import sys

rtod = 180.0/3.141592

def get_args(argv):

    filelist = []
    IsLog = 0
    var = 15
    alt = 400.0
    lon = -100.0
    lat = -100.0
    help = 0


    for arg in argv:

        IsFound = 0

        if (not IsFound):

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = int(m.group(1))
                IsFound = 1


            m = re.match(r'-lat=(.*)',arg)
            if m:
                lat = int(m.group(1))
                IsFound = 1

            m = re.match(r'-lon=(.*)',arg)
            if m:
                lon = int(m.group(1))
                IsFound = 1


            m = re.match(r'-alog',arg)
            if m:
                IsLog = 1
                IsFound = 1

            m = re.match(r'-h',arg)
            if m:
                help = 1
                IsFound = 1

            if IsFound==0 and not(arg==argv[0]):
                filelist.append(arg)

    args = {'filelist':filelist,
            'var':var,
            'help':help,
            'lat':lat,
            'lon':lon,
            'IsLog':IsLog}

    return args

args = get_args(sys.argv)
header = read_gitm_header(args["filelist"])

plon = args['lon']
plat = args['lat']

if (args["help"]):

    print('Usage : ')
    print('gitm_plot_one_alt.py -var=N -lat=lat -lon=lon -alog')
    print('                     -help [file]')
    print('   -help : print this message')
    print('   -var=number : number is variable to plot')
    print('   -lat=latitude : latitude in degrees (closest)')
    print('   -lon=longitude: longitude in degrees (closest)')
    print('   -alog : plot the log of the variable')
    print('   Non-KW args: files.')

    iVar = 0
    for var in header["vars"]:
        print(iVar,var)
        iVar=iVar+1

    exit()

filelist = args["filelist"]
nFiles = len(filelist)
if nFiles < 2:
    print('Please enter multiple files')

iSZA = 43
vars = [0,1,2,iSZA]
vars.append(args["var"])
Var = header["vars"][args["var"]]
AllData2D = []
AllAlts = []
AllTimes = []
AllSZA = []
j = 0
for file in filelist:

    data = read_gitm_one_file(file, vars)
    if (j == 0):
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0]/1000.0;
        Lons = data[0][:,0,0]*rtod;
        Lats = data[1][0,:,0]*rtod;

    ilon = find_nearest_index(Lons,plon)
    ilat = find_nearest_index(Lons,plat)

    AllTimes.append(data["time"])
    AllData2D.append(data[args["var"]][ilon,ilat,:])
    AllSZA.append(data[iSZA][:,:,0])
    AllSZA = np.array(AllSZA)
    pp.figure()
    cb2 = pp.contourf(Lons,Lats,np.transpose(AllSZA[0][:,:]),levels=30)
    imin1 = find_nearest_index(AllSZA[0], 0)
    
    pp.colorbar(cb2)
    pp.savefig('sza.png')
    breakpoint()
AllData2D = np.array(AllData2D) 
AllSZA = np.array(AllSZA)
ialt1 = find_nearest_index(Alts,90)
ialt2 = find_nearest_index(Alts,300)

Alts = np.array(Alts[ialt1:ialt2+1])
AllData2D = AllData2D[:,ialt1:ialt2+1]
fig, ax = plt.subplots()
cont1 = pp.contourf(AllTimes,Alts,np.transpose(AllData2D),levels=30,cmap='turbo')
cb1 = pp.colorbar(cont1,label="{}".format(Var))
myFmt = mdates.DateFormatter("%b %d %H:%M")
ax.xaxis.set_major_formatter(myFmt)

pp.xlabel('Time (UT)')
pp.ylabel('Altitude')
pp.ylim([90,300])
fig.autofmt_xdate()
pp.savefig('plot.png')