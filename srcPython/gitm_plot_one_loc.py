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
    diff = '0'
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

            m = re.match(r'-diff=(.*)',arg)
            if m:
                diff = m.group(1)
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
            'smax':smax,
            'diff':diff}

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
    print('   -alog: plot the log of the variable')
    print('   -diff=backgroundFiles: plot the difference between 2 sets of files')
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
    exit(1)
try:
    iSZA = header["vars"].index('SolarZenithAngle')
    vars = [0,1,2,iSZA]
except: 
    vars = [0,1,2]

diff = False
if args['diff'] != '0':
    diff = True
    backgroundFilelist = sorted(glob(args["diff"]))
    nBackFiles = len(backgroundFilelist)
    if nBackFiles != nFiles:
        print('Difference between sizes of perturbation and background filelists:')
        print('Lengths: {}   {}'.format(nFiles,nBackFiles))
        exit(1)


vars.extend([int(v) for v in args["var"].split(',')])
Var = [header['vars'][int(i)] for i in args['var'].split(',')]
nvars = len(args['var'].split(','))

#We want to store data for multiple variables, so we use a dict where var indices are the keys
AllData = {a:[] for a in args['var'].split(',')}
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
        ialt2 = find_nearest_index(Alts,250)

    AllTimes.append(data["time"])

    if diff:
        stime = str(AllTimes[-1].year)[2:]+str(AllTimes[-1].month).rjust(2,'0')+str(AllTimes[-1].day).rjust(2,'0')+\
            '_'+str(AllTimes[-1].hour).rjust(2,'0')+str(AllTimes[-1].minute).rjust(2,'0')
        bFile = [i for i in backgroundFilelist if stime in i][0]
        if bFile == '':
            #It is possible that we don't have an output file at the same time.
            print('Missing background file corresponding to: {}'.format(file))
            exit(1)
        background = read_gitm_one_file(bFile,vars)
        
    if args['cut'] == 'loc':
        ilon = find_nearest_index(Lons,plon)
        ilat = find_nearest_index(Lons,plat)
        for ivar in args['var'].split(','):
            if diff:
                temp = (data[int(ivar)][ilon,ilat,ialt1:ialt2+1]-background[int(ivar)][ilon,ilat,ialt1:ialt2+1])/ \
                    background[int(ivar)][ilon,ilat,ialt1:ialt2+1]*100.0
            else:
                temp = data[int(ivar)][ilon,ilat,ialt1:ialt2+1]          

            AllData[ivar].append(temp)


    if args['cut'] == 'sza':        
        AllSZA.append(data[iSZA][:,:,0])
        mask = (AllSZA[-1] >= smin) & (AllSZA[-1] <= smax ) 
        for ivar in args['var'].split(','):
            if diff:
                #Calculate the mean of both sets of data and then calculate the percent difference.
                mean1 = data[int(ivar)][:,:,ialt1:ialt2+1][mask].mean(axis=0)
                mean2 = background[int(ivar)][:,:,ialt1:ialt2+1][mask].mean(axis=0)
                temp = (mean1-mean2)/mean2*100.

            else:
                temp = data[int(ivar)][:,:,ialt1:ialt2+1][mask].mean(axis=0)

            AllData[ivar].append(temp)
            # AllData[ivar].append(temp[mask].mean(axis=0))
            
    
    
    j+=1



for ivar in args['var'].split(','):
    AllData[ivar] = np.array(AllData[ivar])

if args['cut']  == 'sza':
    AllSZA = np.array(AllSZA)

# fig, ax = pp.subplots()
fig = pp.figure(figsize=(8.5,11))
pp.ylim([90,300])

Alts = Alts[ialt1:ialt2+1]

cmap = 'plasma'
i=0
for ivar in args['var'].split(','):
    ax = pp.subplot(6,1,i+1)
    AllData2D = AllData[ivar]
<<<<<<< HEAD
    
    if ivar == '3':
        AllData2D = np.log10(AllData2D)
        Var[i] = "Log "+ Var[i]
    
    cont = ax.contourf(AllTimes,Alts,np.transpose(AllData2D),levels=30,cmap='turbo')    
    pp.colorbar(cont,ax=ax,label=Var[i])
    
    #Only plot x axis for last plot
=======
    if ivar == '3' and (not diff):
        AllData2D = np.log10(AllData2D)
        Var[i] = "Log "+ Var[i]

    cont = ax.contourf(AllTimes,Alts,np.transpose(AllData2D),levels=30,cmap=cmap)    
    if diff:
        label = '{}\n% Diff'.format(Var[i])
    else:
        label = Var[i]

    pp.colorbar(cont,ax=ax,label=label)
>>>>>>> 4ce44be46dc1a1092ea4274419b523cc03c57810
    if i < len(Var)-1:
        ax.get_xaxis().set_ticklabels([])
    
    pp.ylabel('Alt (km)')

    i+=1

pp.xlabel('Time (UT)')
myFmt = mdates.DateFormatter("%b %d %H:%M")
ax.xaxis.set_major_formatter(myFmt)
fig.autofmt_xdate()

pp.savefig('plot.png')