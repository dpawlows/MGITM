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
    hmf2 = 0


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

            m = re.match(r'-epeak',arg)
            if m:
                hmf2 = 1
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
            'diff':diff,
            'hmf2':hmf2,
            }

    return args

args = get_args(sys.argv)
header = read_gitm_header(args["filelist"])

if args['hmf2']:
    hmf2 = 1
else:
    hmf2 = 0

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
    print('   -epeak: plot hmf2 as a function of time (will overwrite any vars')
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
altprofile = False
if nFiles < 2:
    altprofile = True
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


if hmf2:
    args['var'] = '32'
    vars.append(int(args['var']))
    nvars = 1
    electron = {args['var']:[],'hmf2':[]}
    Var = [header['vars'][int(args['var'])]]

else:
    vars.extend([int(v) for v in args["var"].split(',')])
    Var = [header['vars'][int(i)] for i in args['var'].split(',')]
    nvars = len(args['var'].split(','))

AllData = {a:[] for a in args['var'].split(',')}

#We want to store data for multiple variables, so we use a dict where var indices are the keys
OCO2 = False
if np.isin(4,vars) and np.isin(6,vars):
    OCO2 = True
    AllData['OCO2'] = []

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

        if OCO2:
            if diff:
                temp = ((data[6][ilon,ilat,ialt1:ialt2+1]/data[4][ilon,ilat,ialt1:ialt2+1]) -\
                    (background[6][ilon,ilat,ialt1:ialt2+1]/background[4][ilon,ilat,ialt1:ialt2+1])) /\
                    (background[6][ilon,ilat,ialt1:ialt2+1]/background[4][ilon,ilat,ialt1:ialt2+1])*100
            else: 
                temp = data[6][ilon,ilat,ialt1:ialt2+1]/data[4][ilon,ilat,ialt1:ialt2+1]

            AllData['OCO2'].append(temp)


    if args['cut'] == 'sza':        
        AllSZA.append(np.array(data[iSZA][2:-2,2:-2,0]))
        
        mask = (AllSZA[-1] >= smin) & (AllSZA[-1] <= smax ) 
     
        for ivar in args['var'].split(','):
            if diff:
                #Calculate the mean of both sets of data and then calculate the percent difference.
                mean1 = data[int(ivar)][2:-2,2:-2,ialt1:ialt2+1][mask].mean(axis=0)
                mean2 = background[int(ivar)][2:-2,2:-2,ialt1:ialt2+1][mask].mean(axis=0)
                temp = (mean1-mean2)/mean2*100.
                if hmf2:
                    perturb = np.amax(data[int(ivar)][2:-2,2:-2,ialt1:ialt2+1],axis=2)[mask].mean()
                    perturbpeak = Alts[ialt1+np.argmax(data[int(ivar)][2:-2,2:-2,ialt1:ialt2+1],axis=2)][mask].mean()
                    back = np.amax(background[int(ivar)][2:-2,2:-2,ialt1:ialt2+1],axis=2)[mask].mean()
                    backpeak = Alts[ialt1+np.argmax(background[int(ivar)][2:-2,2:-2,ialt1:ialt2+1],axis=2)][mask].mean()

                    electron[ivar].append((perturb-back)/back*100)
                    electron['hmf2'].append(perturbpeak)

            else:
                temp = data[int(ivar)][2:-2,2:-2,ialt1:ialt2+1][mask].mean(axis=0)

            AllData[ivar].append(temp)

        if OCO2:
            if diff:
                temp = ((data[6][2:-2,2:-2,ialt1:ialt2+1]/data[4][2:-2,2:-2,ialt1:ialt2+1])[mask].mean(axis=0) -\
                    (background[6][2:-2,2:-2,ialt1:ialt2+1]/background[4][2:-2,2:-2,ialt1:ialt2+1])[mask].mean(axis=0)) /\
                    (background[6][2:-2,2:-2,ialt1:ialt2+1]/background[4][2:-2,2:-2,ialt1:ialt2+1])[mask].mean(axis=0)*100
            else: 
                temp = data[6][2:-2,2:-2,ialt1:ialt2+1]/data[4][2:-2,2:-2,ialt1:ialt2+1]

            AllData['OCO2'].append(temp)

    j+=1


if OCO2:
    args['var'] = args['var']+',OCO2'
    Var.append('O/CO$_2$ ratio')

for ivar in args['var'].split(','):
    AllData[ivar] = np.array(AllData[ivar])

if args['cut']  == 'sza':
    AllSZA = np.array(AllSZA)


if not altprofile:
    fig = pp.figure(figsize=(8.5,11))
    pp.ylim([90,300])

    Alts = Alts[ialt1:ialt2+1]
    cmap = 'plasma'
    # cmap = 'YlOrRd'
    # cmap = 'Spectral'
    maxs = [7.5,39.5,83.2,24.6,55.7,283.5,0]
    # mins = [-.7,-3,-,0,0,0]

    i=0

    for ivar in args['var'].split(','):
        ax = pp.subplot(7,1,i+1)
        AllData2D = AllData[ivar]
        if ivar == '3' and (not diff):
            AllData2D = np.log10(AllData2D)
            Var[i] = "Log "+ Var[i]

            
        if ('maxs' not in locals()):
            maxv = np.amax(AllData2D)
        else:
            maxv = maxs[i]        
        minv = np.amin(AllData2D)

        if ivar == 'OCO2':
            cmap = 'bone'
            minv = -38
        levels = np.linspace(minv,maxv,30)

        cont = ax.contourf(AllTimes,Alts,np.transpose(AllData2D),levels=levels,cmap=cmap)   
        
        if 'single' in filelist[0]:
            line = [datetime(2017, 9, 10, 16, 10)]
        
        if 'doubleflare' in filelist[0]:
            line = [datetime(2017, 9, 10, 16, 10),datetime(2017,9,10,18,45)]

        if 'doubleshort' in filelist[0]:
            line = [datetime(2017, 9, 10, 16, 10),datetime(2017,9,10,16,41)]

        if 'triple' in filelist[0]:
            line = [datetime(2017, 9, 10, 16, 10),datetime(2017,9,10,16,41),datetime(2017,9,10,17,11)]

        line2 = datetime(2017, 9, 10, 23, 30)

        if 'line' in locals():
            for l in line:
                ax.plot([l,l],[Alts[0],Alts[-1]],color='dimgrey',linestyle=(0,(10,3)))

            ax.plot([line2,line2],[Alts[0],Alts[-1]],color='dimgrey',linestyle=(0,(5,10)))

        if Var[i] in name_dict:
            varname = name_dict[Var[i]]
        else:
            varname = Var[i]

        label = varname
        if diff:
            label = label+'\n% Diff'

        pp.colorbar(cont,ax=ax,label=label,format="%.1f")
        if i < len(Var)-1:
            ax.get_xaxis().set_ticklabels([])
        pp.ylabel('Alt (km)')
        i+=1

    pp.xlabel('Time (UT)')
    myFmt = mdates.DateFormatter("%b %d %H:%M")
    ax.xaxis.set_major_formatter(myFmt)
    fig.autofmt_xdate()
    pp.savefig('plot.png')

    if hmf2:
        fig1 = pp.figure()
        ax = pp.subplot(2,1,1)
        # ax.get_xaxis().set_ticklabels([])

        ax.plot(AllTimes,electron['32'])
        ax.set_ylim([0,22])
        pp.ylabel('Peak [e-] % Diff')

        ax1 = pp.subplot(2,1,2)
        ax1.plot(AllTimes,electron['hmf2'])
        pp.ylabel('Alt of Peak [e-]')
        pp.xlabel('Time (UT)')
        myFmt = mdates.DateFormatter("%b %d %H:%M")
        ax1.xaxis.set_major_formatter(myFmt)
        fig1.autofmt_xdate()
        pp.savefig('electron.png')

else:
    #Then we are only plotting a single file and it should be an altitude profile
    fig = pp.figure(figsize=(8.5,11))
    pp.ylim([90,300])

    Alts = Alts[ialt1:ialt2+1]
    i = 0
    for ivar in args['var'].split(','):
        ax = pp.subplot(2,6,i+1)
        AllData2D = AllData[ivar]
        if ivar == '3' and (not diff):
            AllData2D = np.log10(AllData2D)
            Var[i] = "Log "+ Var[i]

            
        if ('maxs' not in locals()):
            maxv = np.amax(AllData2D)
        else:
            maxv = maxs[i]        
        minv = np.amin(AllData2D)
        ax.plot(AllData2D[0,:],Alts)
        pp.xlim([0,325])
        if Var[i] in name_dict:
            varname = name_dict[Var[i]]
        else:
            varname = Var[i]

        label = varname
        if diff:
            label = label+'\n% Diff'
        pp.xlable = label
        i+=1


    pp.ylabel('Alt (km)')

    pp.savefig('plot.png')


### One time use plot a single altitude as a function of time (or output data)
# fig1 = pp.figure()
# lineplot = '3'
# palt = 200
# ialt = find_nearest_index(Alts, palt)
# ax2 = pp.subplot(2,1,1)

# if 'lineplot' in locals():
#     ax2.plot(AllTimes,AllData[lineplot][:,ialt])

# pp.ylim([0,33])
# pp.xlabel('Time (UT)')
# pp.ylabel('{}\n% Diff ({} km)'.format(name_dict[Var[vars.index(int(lineplot))-4]],int(Alts[ialt])))
# myFmt = mdates.DateFormatter("%b %d %H:%M")
# ax2.xaxis.set_major_formatter(myFmt)
# fig1.autofmt_xdate()
# pp.savefig('lineplot.png')

filename = 'output.txt'
f = open(filename,'a')
o1 = [a.strftime("%d-%b-%y/%H:%M:%S") for a in AllTimes]
output = o1  + list(electron['32'])
for o in output:
    f.write(str(o)+' ')
f.write("\n")
f.close()

