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
    cut = 'loc'
    smin = -100.0
    smax = -100.0
    help = 0


    for arg in argv:

        IsFound = 0

        if (not IsFound):

            m = re.match(r'-var=(.*)',arg)
            if m:
                var = m.group(1)
                IsFound = 1

            m = re.match(r'-cut=(.*)',arg)
            if m:
                cut = m.group(1)
                IsFound = 1


            m = re.match(r'-lat=(.*)',arg)
            if m:
                lat = int(m.group(1))
                IsFound = 1

            m = re.match(r'-lon=(.*)',arg)
            if m:
                lon = int(m.group(1))
                IsFound = 1

            m = re.match(r'-smin=(.*)',arg)
            if m:
                smin = int(m.group(1))
                IsFound = 1   

            m = re.match(r'-smax=(.*)',arg)
            if m:
                smax = int(m.group(1))
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
            'IsLog':IsLog,
            'cut':cut,
            'smin':smin,
            'smax':smax,}

    return args

args = get_args(sys.argv)
header = read_gitm_header(args["filelist"])

if args['cut'] == 'loc' and args['lon'] > -50:
    plon = args['lon']
    plat = args['lat']
elif args['cut'] == 'sza' and args['smin'] > -50:
    smin = args['smin']
    smax = args['smax']
else:
    args["help"] = '-h'


if (args["help"]):

    print('Usage : ')
    print('gitm_plot_one_alt.py -var=N1[,N2,N3,...] -lat=lat -lon=lon -alog')
    print('                     -help [file]')
    print('   -help : print this message')
    print('   -var=number[,num2,num3,...] : number is variable to plot')
    print('   -cut=loc,sza: Plot type ')
    print('   -lat=latitude : latitude in degrees (closest) (cut=loc) ')
    print('   -lon=longitude: longitude in degrees (closest) (cut=loc)')
    print('   -smin=minsza: minimum solar zenith angle (cut=sza)')
    print('   -smax=maxsza: maximum solar zenigh angle (cut=sza)')
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
try:
    iSZA = header["vars"].index('SolarZenithAngle')
    vars = [0,1,2,iSZA]
except: 
    vars = [0,1,2]

vars.extend([int(v) for v in args["var"].split(',')])
Var = [header['vars'][int(i)] for i in args['var'].split(',')]

#We want to store data for multiple variables, so we use a dict where var indices are the keys
AllData = {a for a in args['var'].split(',')}
AllData2D = []
AllAlts = []
AllTimes = []
AllSZA = []
j = 0
for file in filelist:

    data = read_gitm_one_file(file, vars)
    if (j == 0):    
        [nLons, nLats, nAlts] = data[0].shape
        Alts = data[2][0][0]/1000.0
        Lons = data[0][:,0,0]*rtod
        Lats = data[1][0,:,0]*rtod

        ialt1 = find_nearest_index(Alts,90)
        ialt2 = find_nearest_index(Alts,300)


    AllTimes.append(data["time"])
    if args['cut'] == 'loc':
        ilon = find_nearest_index(Lons,plon)
        ilat = find_nearest_index(Lons,plat)
        AllData2D.append(data[args["var"]][ilon,ilat,:])

    if args['cut'] == 'sza':        
        AllSZA.append(data[iSZA][:,:,0])
        breakpoint()
        for int(ivar) in args['var'].split(','):
            temp = 
        mask = (AllSZA[-1] >= smin) & (AllSZA[-1] <= smax ) 
        temp = data[args["var"]][:,:,ialt1:]
        AllData2D.append(temp[mask].mean(axis=0))
     
        j+=1

AllData2D = np.array(AllData2D) 
fig, ax = pp.subplots()
myFmt = mdates.DateFormatter("%b %d %H:%M")
ax.xaxis.set_major_formatter(myFmt)


if args['cut']  == 'loc':
    Alts = Alts[ialt1:ialt2+1]
    AllData2D = AllData2D[:,ialt1:ialt2+1]
    cont1 = pp.contourf(AllTimes,Alts,np.transpose(AllData2D),levels=30,cmap='turbo')
    cb1 = pp.colorbar(cont1,label="{}".format(Var))


if args['cut']  == 'sza':
    AllSZA = np.array(AllSZA)
    Alts = Alts[ialt1:]

    cont1 = pp.contourf(AllTimes,Alts,np.transpose(AllData2D),levels=30,cmap='turbo')
    cb1 = pp.colorbar(cont1,label="{}".format(Var))


pp.xlabel('Time (UT)')
pp.ylabel('Altitude')
pp.ylim([90,300])
fig.autofmt_xdate()


pp.savefig('plot.png')