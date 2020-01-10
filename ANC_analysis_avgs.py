# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:46:15 2020

@author: admin
"""

import pandas as pd
import numpy as np

import dill

import matplotlib.pyplot as plt

def calcavg(df, day, col):
    
    """ Calculates mean and SEM from specified columns of dataframe for phase 1 and phase 2
    
    Parameters:
        df: dataframe with data from all days for individual rats in ANC expt
        day: day to return data for - MUST be a string!
        col: column to calculate without phase suffix
        
    Returns:
        x: list of mean values for phase 1 and phase 2
        sem: list of sem values for phase 1 and phase 2
    """
    
    df = df.loc[df['day'] == day]
    
    x = []
    sem = []
    
    for ph in ['_ph1', '_ph2']:
        x.append(np.mean(df[col+ph]))
        sem.append(np.std(df[col+ph]) / len(df[col+ph]))
    
    return x, sem


folder='C:\\Github\\anc-for-jas\\data\\jess\\'
        
pickle_in = open(folder + 'anc_dataframe.pickle', 'rb')
df = dill.load(pickle_in)

df_avg = pd.DataFrame()

days = df.day.unique()

#for day in days:
#    df_avg['mean']


"""
Next step is to use matplotlib to make figure
"""

f, ax = plt.subplots(nrows=7, ncols=2, figsize=(10,12), sharey=True, sharex=True)

figcol=0
for idx, day in enumerate(['1.0', '3.0', '5.0', '7.0', '9.0', '11.0', '13.0']):
    mean, sem = calcavg(df, day, 'nLicks')
    axis=ax[idx][figcol]
    axis.bar([1,2], mean)
    axis.set_ylim([-20, 2000])
    axis.set_xticks([1, 2])
    axis.set_xticklabels(['Phase 1', 'Phase 2'])
    axis.set_ylabel('Number of Licks')

    
figcol=1
for idx, day in enumerate(['2.0', '4.0', '6.0', '8.0', '10.0', '12.0', '14.0']):
    mean, sem = calcavg(df, day, 'nLicks')
    axis=ax[idx][figcol]
    axis.bar([1,2], mean, color='orange')
    axis.set_ylim([-20, 2000])
    
f.savefig('C:\\Github\\anc-for-jas\\data\\jess\\figs\\quick nLicks figure.png')


