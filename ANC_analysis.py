# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:33:13 2019

@author: admin
"""

# ANC analysis for Celia and John

# import packages
import numpy as np
import scipy.optimize as opt
import scipy.stats as stats

import matplotlib.pyplot as plt

# Necessary functions
def medfilereader(filename, varsToExtract = 'all',
                  sessionToExtract = 1,
                  verbose = False,
                  remove_var_header = False):
    if varsToExtract == 'all':
        numVarsToExtract = np.arange(0,26)
    else:
        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
    
    f = open(filename, 'r')
    f.seek(0)
    filerows = f.readlines()[8:]
    datarows = [isnumeric(x) for x in filerows]
    matches = [i for i,x in enumerate(datarows) if x == 0.3]
    if sessionToExtract > len(matches):
        print('Session ' + str(sessionToExtract) + ' does not exist.')
    if verbose == True:
        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
        print('Analyzing session ' + str(sessionToExtract))
    
    varstart = matches[sessionToExtract - 1]
    medvars = [[] for n in range(26)]
    
    k = int(varstart + 27)
    for i in range(26):
        medvarsN = int(datarows[varstart + i + 1])
        
        medvars[i] = datarows[k:k + int(medvarsN)]
        k = k + medvarsN
        
    if remove_var_header == True:
        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
    else:
        varsToReturn = [medvars[i] for i in numVarsToExtract]

    if np.shape(varsToReturn)[0] == 1:
        varsToReturn = varsToReturn[0]
    return varsToReturn


def isnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')

def lickCalc(licks, offset = [], burstThreshold = 0.5, runThreshold = 10,
             ignorelongilis=True, minburstlength=1,
             binsize=60, histDensity = False,
             perform_Weibull=False):
    # makes dictionary of data relating to licks and bursts
    if type(licks) != np.ndarray or type(offset) != np.ndarray:
        try:
            licks = np.array(licks)
            offset = np.array(offset)
        except:
            print('Licks and offsets need to be arrays and unable to easily convert.')
            return

    lickData = {}
    
    if len(offset) > 0:
        lickData['licklength'] = offset - licks
        lickData['longlicks'] = [x for x in lickData['licklength'] if x > 0.3]
    else:
        lickData['licklength'] = []
        lickData['longlicks'] = []

    lickData['licks'] = np.concatenate([[0], licks])
    lickData['ilis'] = np.diff(lickData['licks'])
    if ignorelongilis:
        lickData['ilis'] = [x for x in lickData['ilis'] if x < burstThreshold]

    lickData['freq'] = 1/np.mean([x for x in lickData['ilis'] if x < burstThreshold])
    lickData['total'] = len(licks)
    
    # Calculates start, end, number of licks and time for each BURST 
    lickData['bStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]  
    lickData['bInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]
    lickData['bEnd'] = [lickData['licks'][i-1] for i in lickData['bInd'][1:]]
    lickData['bEnd'].append(lickData['licks'][-1])

    lickData['bLicks'] = np.diff(lickData['bInd'] + [len(lickData['licks'])])
    
    # calculates which bursts are above the minimum threshold
    inds = [i for i, val in enumerate(lickData['bLicks']) if val>minburstlength-1]
    
    lickData['bLicks'] = removeshortbursts(lickData['bLicks'], inds)
    lickData['bStart'] = removeshortbursts(lickData['bStart'], inds)
    lickData['bEnd'] = removeshortbursts(lickData['bEnd'], inds)
      
    lickData['bTime'] = np.subtract(lickData['bEnd'], lickData['bStart'])
    lickData['bNum'] = len(lickData['bStart'])
    if lickData['bNum'] > 0:
        lickData['bMean'] = np.nanmean(lickData['bLicks'])
        lickData['bMean-first3'] = np.nanmean(lickData['bLicks'][:3])
    else:
        lickData['bMean'] = 0
        lickData['bMean-first3'] = 0
    
    lickData['bILIs'] = [x for x in lickData['ilis'] if x > burstThreshold]

    # Calculates start, end, number of licks and time for each RUN
    lickData['rStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]  
    lickData['rInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]
    lickData['rEnd'] = [lickData['licks'][i-1] for i in lickData['rInd'][1:]]
    lickData['rEnd'].append(lickData['licks'][-1])

    lickData['rLicks'] = np.diff(lickData['rInd'] + [len(lickData['licks'])])    
    lickData['rTime'] = np.subtract(lickData['rEnd'], lickData['rStart'])
    lickData['rNum'] = len(lickData['rStart'])

    lickData['rILIs'] = [x for x in lickData['ilis'] if x > runThreshold]
    if perform_Weibull:
        try:
            xdata, ydata = calculate_burst_prob(lickData['bLicks'])
            try:
                lickData['weib_alpha'], lickData['weib_beta'], lickData['weib_rsq'] = fit_weibull(xdata, ydata)
            except RuntimeError:
                print('Optimal fit parameters not found')
                lickData['weib_alpha'], lickData['weib_beta'], lickData['weib_rsq'] = [0, 0, 0]
                
            lickData['burstprob']=[xdata, ydata]
        except ValueError:
            print('Could not calculate burst probability')
            lickData['weib_alpha'], lickData['weib_beta'], lickData['weib_rsq'] = [0, 0, 0]
            lickData['burstprob']=[[0], [0]]
    try:
        lickData['hist'] = np.histogram(lickData['licks'][1:], 
                                    range=(0, 3600), bins=int((3600/binsize)),
                                          density=histDensity)[0]
    except TypeError:
        print('Problem making histograms of lick data')
        
    return lickData

def removeshortbursts(data, inds):
    data = [data[i] for i in inds]
    return data

def calculate_burst_prob(bursts):
    bins = np.arange(min(bursts), max(bursts))
    hist=np.histogram(bursts, bins=bins, density=True)
    cumsum=np.cumsum(hist[0])

    x = hist[1][1:]
    y = [1-val for val in cumsum]
    
    return x, y

def weib_davis(x, alpha, beta): 
    return (np.exp(-(alpha*x)**beta))

def fit_weibull(xdata, ydata):
    x0=np.array([0.1, 1])
    fit=opt.curve_fit(weib_davis, xdata, ydata, x0)
    alpha=fit[0][0]
    beta=fit[0][1]
    slope, intercept, r_value, p_value, std_err = stats.linregress(ydata, weib_davis(xdata, alpha, beta))
    r_squared=r_value**2
    
    return alpha, beta, r_squared

def barscatter(data, transpose = False, unequal=False,
                groupwidth = .75,
                barwidth = .9,
                paired = False,
                spaced = False,
                yspace = 20,
                xspace = 0.1,
                barfacecoloroption = 'same', # other options 'between' or 'individual'
                barfacecolor = ['white'],
                baredgecoloroption = 'same',
                baredgecolor = ['black'],
                baralpha = 1,
                scatterfacecoloroption = 'same',
                scatterfacecolor = ['white'],
                scatteredgecoloroption = 'same',
                scatteredgecolor = ['grey'],
                scatterlinecolor = 'grey', # Don't put this value in a list
                scattersize = 80,
                scatteralpha = 1,
                spreadscatters = False,
                linewidth=1,
                xlim=[],
                ylim=[],
                ylabel = 'none',
                xlabel = 'none',
                grouplabel = 'auto',
                itemlabel = 'none',
                barlabels = [],
                barlabeloffset=0.025,
                grouplabeloffset=0,
                yaxisparams = 'auto',
                show_legend = 'none',
                legendloc='upper right',
                xfontsize=8,
                ax=[]):

    if unequal == True:
        dims = np.ndim(data)
        data_obj = np.ndarray((np.shape(data)), dtype=np.object)
        if dims == 1:
            for i, dim in enumerate(data):
                data_obj[i] = np.array(dim, dtype=np.object)
            data = data_obj
        elif dims == 2:            
            for i1, dim1 in enumerate(data):
                for i2, dim2 in enumerate(dim1):
                    data_obj[i1][i2] = np.array(dim2, dtype=np.object)
            data = data_obj
        else:
            print('Cannot convert that number of dimensions or data is in wrong format. Attmepting to make graph assuming equal groups.')
    
    if type(data) != np.ndarray or data.dtype != np.object:
        dims = np.shape(data)
        if len(dims) == 2 or len(dims) == 1:
            data = data2obj1D(data)

        elif len(dims) == 3:
            data = data2obj2D(data)
              
        else:
            print('Cannot interpret data shape. Should be 2 or 3 dimensional array. Exiting function.')
            return

    # Check if transpose = True
    if transpose == True:
        data = np.transpose(data)
        
    # Initialize arrays and calculate number of groups, bars, items, and means
    
    barMeans = np.zeros((np.shape(data)))
    items = np.zeros((np.shape(data)))
    
    nGroups = np.shape(data)[0]
    groupx = np.arange(1,nGroups+1)

    if len(np.shape(data)) > 1:
        grouped = True
        barspergroup = np.shape(data)[1]
        barwidth = (barwidth * groupwidth) / barspergroup
        
        for i in range(np.shape(data)[0]):
            for j in range(np.shape(data)[1]):
                barMeans[i][j] = np.mean(data[i][j])
                items[i][j] = len(data[i][j])
        
    else:
        grouped = False
        barspergroup = 1
        
        for i in range(np.shape(data)[0]):
            barMeans[i] = np.mean(data[i])
            items[i] = len(data[i])
    
    # Calculate x values for bars and scatters
    
    xvals = np.zeros((np.shape(data)))
    barallocation = groupwidth / barspergroup
    k = (groupwidth/2) - (barallocation/2)
    
    if grouped == True:
        
        for i in range(np.shape(data)[0]):
            xrange = np.linspace(i+1-k, i+1+k, barspergroup)
            for j in range(barspergroup):
                xvals[i][j] = xrange[j]
    else:
        xvals = groupx
    
    # Set colors for bars and scatters
     
    barfacecolorArray = setcolors(barfacecoloroption, barfacecolor, barspergroup, nGroups, data)
    baredgecolorArray = setcolors(baredgecoloroption, baredgecolor, barspergroup, nGroups, data)
     
    scfacecolorArray = setcolors(scatterfacecoloroption, scatterfacecolor, barspergroup, nGroups, data, paired_scatter = paired)
    scedgecolorArray = setcolors(scatteredgecoloroption, scatteredgecolor, barspergroup, nGroups, data, paired_scatter = paired)
    
    # Initialize figure
    if ax == []:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    # Make bars
    barlist = []
    barx = []
    for x, y, bfc, bec in zip(xvals.flatten(), barMeans.flatten(),
                              barfacecolorArray, baredgecolorArray):
        barx.append(x)
        barlist.append(ax.bar(x, y, barwidth,
                         facecolor = bfc, edgecolor = bec,
                         zorder=-1))
    
    # Uncomment these lines to show method for changing bar colors outside of
    # function using barlist properties
    #for i in barlist[2].get_children():
    #    i.set_color('r')
    
    # Make scatters
    sclist = []
    if paired == False:
        for x, Yarray, scf, sce  in zip(xvals.flatten(), data.flatten(),
                                        scfacecolorArray, scedgecolorArray):
            if spaced == True:
                xVals, yVals = xyspacer(ax, x, Yarray, bindist=yspace, space=xspace)
                sclist.append(ax.scatter(xVals, yVals, s = scattersize,
                             c = scf,
                             edgecolors = sce,
                             zorder=20))           
            else:
                for y in Yarray:
                     sclist.append(ax.scatter(x, y, s = scattersize,
                                     c = scf,
                                     edgecolors = sce,
                                     zorder=20))

    elif grouped == True:
        for x, Yarray, scf, sce in zip(xvals, data, scfacecolorArray, scedgecolorArray):
            for y in np.transpose(Yarray.tolist()):
                sclist.append(ax.plot(x, y, '-o', markersize = scattersize/10,
                         color = scatterlinecolor,
                         linewidth=linewidth,
                         markerfacecolor = scf,
                         markeredgecolor = sce,
                         zorder=20))
    elif grouped == False:
        for n,_ in enumerate(data[0]):
            y = [y[n-1] for y in data]
            sclist.append(ax.plot(xvals, y, '-o', markersize = scattersize/10,
                         color = scatterlinecolor,
                         linewidth=linewidth,
                         markerfacecolor = scfacecolorArray[0],
                         markeredgecolor = scedgecolorArray[0],
                         zorder=20))
    
    # Label axes
    if ylabel != 'none':
        ax.set_ylabel(ylabel)
    
    if xlabel != 'none':
        ax.set_xlabel(xlabel)
    
    # Set range and tick values for Y axis
    if yaxisparams != 'auto':
        ax.set_ylim(yaxisparams[0])
        ax.set_yticks(yaxisparams[1])
       
    # X ticks
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off') # labels along the bottom edge are off

    if len(xlim) > 0:
        ax.set_xlim(xlim)
    if len(ylim) > 0:
        ax.set_ylim(ylim)
        
    xrange = ax.get_xlim()[1] - ax.get_xlim()[0]
    yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
        
    if grouplabel == 'auto':
        ax.tick_params(labelbottom='off')
    else:
        ax.tick_params(labelbottom='off')

        groupx = np.arange(1, len(grouplabel)+1)
        if len(xlim) > 0:
            groupx = [x for x in groupx]
        xpos = (groupx - ax.get_xlim()[0])/xrange

        for x, label in zip(xpos, grouplabel):
            ax.text(x, -0.05+grouplabeloffset, label, va='top', ha='center', fontsize=xfontsize, transform=ax.transAxes)
        
    if len(barlabels) > 0:
        if len(barlabels) != len(barx):
            print('Wrong number of bar labels for number of bars!')
        else:
            xpos = (barx - ax.get_xlim()[0])/xrange
            ypos = (-ax.get_ylim()[0]/yrange) - barlabeloffset

            for x, label in zip(xpos, barlabels):
                ax.text(x, ypos, label, va='top', ha='center', fontsize=xfontsize, transform=ax.transAxes)

    # Hide the right and top spines and set bottom to zero
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    
    if show_legend == 'within':
        if len(itemlabel) != barspergroup:
            print('Not enough item labels for legend!')
        else:
            legendbar = []
            legendtext = []
            for i in range(barspergroup):
                legendbar.append(barlist[i])
                legendtext.append(itemlabel[i])
            plt.legend(legendbar, legendtext, loc=legendloc)
    
    return ax, barx, barlist, sclist

#plt.savefig('foo.png')
        
# To do
# check if n's are the same for paired and if not default to unpaired
# add color options for scatters
# add alpha options etc
# add axis options
# remove white background
# work out how to export or save as pdf, tiff, eps etc
# work out how to return handles to scatters so can be altered outside of function
# make help doc
# make html file to show usage using ijupyter

      
def setcolors(coloroption, colors, barspergroup, nGroups, data, paired_scatter = False):
            
    nColors = len(colors)
    
    if (paired_scatter == True) & (coloroption == 'within'):
        print('Not possible to make a Paired scatter plot with Within setting.')
        coloroption = 'same'
        
    if coloroption == 'within':
        if nColors < barspergroup:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > barspergroup:
            colors = colors[:barspergroup]
        coloroutput = [colors for i in data]
        coloroutput = list(chain(*coloroutput))
        
    if coloroption == 'between':
        if nColors < nGroups:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > nGroups:
            colors = colors[:nGroups]
        if paired_scatter == False:
            coloroutput = [[c]*barspergroup for c in colors]
            coloroutput = list(chain(*coloroutput))
        else:
            coloroutput = colors
            
    if coloroption == 'individual':
        if nColors < nGroups*barspergroup:
            print('Not enough colors for this color option')
            coloroption = 'same'
        elif nColors > nGroups*barspergroup:
            coloroutput = colors[:nGroups*barspergroup]
        else: 
            coloroutput = colors
    
    if coloroption == 'same':
        coloroutput = [colors[0] for x in range(len(data.flatten()))]

    return coloroutput

def data2obj1D(data):
    obj = np.empty(len(data), dtype=np.object)
    for i,x in enumerate(data):
        obj[i] = np.array(x)  
    return obj

def data2obj2D(data):
    obj = np.empty((np.shape(data)[0], np.shape(data)[1]), dtype=np.object)
    for i,x in enumerate(data):
        for j,y in enumerate(x):
            obj[i][j] = np.array(y)
    return obj

def shadedError(ax, yarray, linecolor='black', errorcolor = 'xkcd:silver', linewidth=1):
    yarray = np.array(yarray)
    y = np.mean(yarray, axis=0)
    yerror = np.std(yarray)/np.sqrt(len(yarray))
    x = np.arange(0, len(y))
    ax.plot(x, y, color=linecolor, linewidth=1)
    ax.fill_between(x, y-yerror, y+yerror, color=errorcolor, alpha=0.4)
    
    return ax

def openANCfile(filename, session):
    
    ph1_licks = {}
    ph1_licks['onset'], ph1_licks['offset'] = medfilereader(filename,
                         varsToExtract=['b', 'c'],
                         sessionToExtract=session,
                         remove_var_header = True)
    
    ph1 = lickCalc(ph1_licks['onset'], burstThreshold = 1, runThreshold = 10,
                        ignorelongilis=True, minburstlength=3,
                        binsize=60, histDensity = False)
    
    
    ph2_licks = {}
    ph2_licks['onset'], ph2_licks['offset'] = medfilereader(filename,
                         varsToExtract=['e', 'f'],
                         sessionToExtract=1,
                         remove_var_header = True)
    
    ph2 = lickCalc(ph2_licks['onset'], burstThreshold = 1, runThreshold = 10,
                         ignorelongilis=True, minburstlength=3,
                         binsize=60, histDensity = False)
    
    h1 = ph1['hist'][:5]
    h2 = ph2['hist'][5:10]   
    hist = np.concatenate((h1, h2), axis=0)

    return ph1, ph2, hist

# Open files

folder = 'C:\\Github\\anc-for-jas\\data\\'

filenames = ['!2019-10-07_09h43m.Subject ANC_day30_C10_LR_1',
             '!2019-10-07_09h51m.Subject ANC_day30_C10_LR_3',
             '!2019-10-07_10h03m.Subject ANC_day30_C10_RL_4',
             '!2019-10-07_10h13m.Subject ANC_day30_C10_RL_6',
             '!2019-10-07_10h38m.Subject ANC_day30_C10_RL_5',
             '!2019-10-07_10h48m.Subject ANC_day30_C10_RL_8',
             '!2019-10-07_11h02m.Subject ANC_day30_C10_LR_2',
             '!2019-10-07_11h12m.Subject ANC_day30_C10_LR_7']

#filenames=['!2019-10-07_10h03m.Subject ANC_day30_C10_RL_4']

Ph1 = {}
Ph2 = {}

for Ph in [Ph1, Ph2]:
    Ph['licks'] = []
    Ph['total'] = []
    Ph['bMean'] = []
    Ph['bNum'] = []
    
anc_hist=[]

for file in filenames:
    for s in [1, 2]:
        ph1, ph2, hist = openANCfile(folder+file, s)
        for Ph, ph in zip([Ph1, Ph2], [ph1, ph2]):
            Ph['licks'].append(ph['licks'])
            Ph['total'].append(ph['total'])
            Ph['bMean'].append(ph['bMean'])
            Ph['bNum'].append(ph['bNum'])
        anc_hist.append(hist)

    
# Set up figure        
f, ax = plt.subplots(figsize=(8, 6), ncols=2, nrows=2)        
f.subplots_adjust(wspace=0.4, hspace=0.4)

# Prepare data for panel A
data = [Ph1['total'], Ph2['total']]

barscatter(data, paired=True, ax=ax[0][0])

ax[0][0].set_ylabel('Total licks')
ax[0][0].set_ylim([-50, 1000])

# Prepare data for panel B

#ax[0][1].plot(np.mean(anc_hist, axis=0))
#ax[0][1].plot(np.transpose(anc_hist))

shadedError(ax[0][1], anc_hist)

ax[0][1].set_ylabel('Licks per min')
ax[0][1].set_xlabel('Time (min)')



data = [Ph1['bNum'], Ph2['bNum']]

barscatter(data, paired=True, ax=ax[1][0])

ax[1][0].set_ylabel('Number of bursts')
ax[1][0].set_ylim([-1, 15])


data = [Ph1['bMean'], Ph2['bMean']]

barscatter(data, paired=True, ax=ax[1][1])
ax[1][1].set_ylabel('Mean licks/burst')
ax[1][1].set_ylim([-50, 400])
        
#        f, ax = plt.subplots(ncols=2)
#        f.subplots_adjust(wspace=0.4)
#        ax[0].bar(['Ph1', 'Ph2'], [ph1['total'], ph2['total']])
#        ax[0].set_ylabel('Total licks')
#        
#        ax[1].bar(['Ph1', 'Ph2'], [ph1['bMean'], ph2['bMean']])
#        ax[1].set_ylabel('Mean licks per burst')






