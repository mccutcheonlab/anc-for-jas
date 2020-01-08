# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:29:49 2020

@author: admin
"""

import pandas as pd
import xlrd

from ANC_helperfx import *

def extract_data(filename, burstThreshold=0.5, minburstlength=1):
    """ Extract data from MedAssociates files in ANC experiment
    
    Parameters:
        filename: Should be complete filename including folder
        burstThreshold: shortest ILI that separates runs of licks
        minburstlength: minimum number of licks that defines a burst
        
    Returns:
        Lick data for left and right bottle
    """

    print('Extracting data from {}'.format(filename))
    
    
    L_licks={}
    R_licks={}
    
    for sidedict, medvars in zip([L_licks, R_licks],
                                 [['b', 'c'], ['e', 'f']]):    
        
        sidedict['onset'], sidedict['offset'] = medfilereader(filename,
                             varsToExtract=medvars,
                             remove_var_header = True)
        
        sidedict['lickdata'] = lickCalc(sidedict['onset'], offset=sidedict['offset'],
                     burstThreshold = burstThreshold,
                     runThreshold = 10,
                     ignorelongilis=True,
                     minburstlength=minburstlength,
                     binsize=60, histDensity = False)
    
    return [L_licks, R_licks]

def assemble_dataframe(lickdata, cols):
    print('To assemble dataframe to be appended to master df')

cols=['filename',
      'nLicks(L)',
      't of first lick(L)',
      '# of long licks',
      'licks per burst',
      '# of bursts']

df = pd.DataFrame(columns=cols)

folder='C:\\Github\\anc-for-jas\\data\\'
filename=folder+'!2019-11-07_15h10m.Subject 0'

data = extract_data(filename, minburstlength=3)

#
#df['filename'] = [filename]
#
#
#filename=folder+'!2019-11-07_15h10m.Subject 1'
#
#df['filename'] = [filename]

#ph1_licks = {}
#ph1_licks['onset'], ph1_licks['offset'] = medfilereader(filename,
#                     varsToExtract=['b', 'c'],
#                     sessionToExtract=session,
#                     remove_var_header = True)
#
#ph1 = lickCalc(ph1_licks['onset'], burstThreshold = 1, runThreshold = 10,
#                    ignorelongilis=True, minburstlength=3,
#                    binsize=60, histDensity = False)