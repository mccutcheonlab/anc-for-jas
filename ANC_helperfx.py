# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:40:42 2020

@author: admin
"""

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
