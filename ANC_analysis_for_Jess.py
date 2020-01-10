# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:29:49 2020

@author: admin
"""
import numpy as np
import pandas as pd
import xlrd as xl
import os

import dill

import matplotlib.pyplot as plt

from ANC_helperfx import *

def extract_data(filename, session=1, burstThreshold=0.5, minburstlength=1):
    """ Extract data from MedAssociates files in ANC experiment
    
    Parameters:
        filename: Should be complete filename including folder
        session: Session to extract (if >1)
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
                             sessionToExtract=session,
                             remove_var_header = True)
        
        sidedict['onset'], sidedict['offset'] = equal_length_lists(sidedict['onset'], sidedict['offset'])
        
        sidedict['lickdata'] = lickCalc(sidedict['onset'], offset=sidedict['offset'],
                     burstThreshold = burstThreshold,
                     runThreshold = 10,
                     ignorelongilis=True,
                     minburstlength=minburstlength,
                     binsize=30, histDensity = False)
    
    return [L_licks, R_licks]

def equal_length_lists(listA, listB):
    if len(listA) == len(listB):
        return listA, listB
    else:
        if len(listA) > len(listB):
            listA = listA[:-1]
        else:
            listB = listB[1:]
        return listA, listB
    

def assemble_dataframe(lickdata, suffix=''):
    """ Assembles lickdata into a dataframe with specified suffix as column name
    
    Parameters:
        lickdata: dictionary of lickdata produced by lickcalc function
        suffix: string specifying suffix to be used when naming columns
        
    Returns:
        Dataframe containign specified values
        
    """
    sf=suffix
    
    df = pd.DataFrame()
    df['nLicks' + sf] = [lickdata['total']]
    if lickdata['total']>0:
        df['t of first lick' + sf] = [lickdata['licks'][1]]
        df['t of last lick' + sf] = [lickdata['licks'][lickdata['total']]]
        df['# of longlicks' + sf] = [len(lickdata['longlicks'])]
        df['licks per burst' + sf] = [lickdata['bMean']]
        df['# of bursts' + sf] = [lickdata['bNum']]
        
    else:
        df['t of first lick' + sf] = [np.nan]
        df['t of last lick' + sf] = [np.nan]
        df['# of longlicks' + sf] = [np.nan]
        df['licks per burst' + sf] = [np.nan]
        df['# of bursts' + sf] = [np.nan]

    return df

def make_histogram_fig(lickdataL, lickdataR, file=''):
    f, ax = plt.subplots()
    ax.plot(lickdataL['hist'], 'blue')
    ax.plot(lickdataR['hist'], 'red')
    
    ax.set_ylabel('Licks per 30s bin')
    
    ax.set_xlim([0, 20])
    ax.set_xticks([5, 15])
    ax.set_xticklabels(['Phase 1', 'Phase 2'])
    
    ax.plot([10, 10], ax.get_ylim(), '--', color='xkcd:charcoal')
    
    ax.set_title(file)
    
    return f

analyseallfilesinfolder=False

folder='C:\\Github\\anc-for-jas\\data\\jess\\'

list_of_files = os.listdir(folder)

df=pd.DataFrame()

### To analyse all files in folder

if analyseallfilesinfolder:

    for file in list_of_files:
        for session in [1,2]:
            try:
                filename=folder+file
                
                data = extract_data(filename, session=session, minburstlength=3)
                
                fig = make_histogram_fig(data[0]['lickdata'], data[1]['lickdata'], file=file+'_session '+str(session))
                fig.savefig('C:\\Github\\anc-for-jas\\data\\jess\\figs\\' + file + '.png')
                
                dfL = assemble_dataframe(data[0]['lickdata'], suffix='_L')
                dfR = assemble_dataframe(data[1]['lickdata'], suffix='_R')
                
                df_id = pd.DataFrame()
                df_id['filename'] = [file]
                df_id['session'] = [str(session)]
                
                df_temp = pd.concat([df_id, dfL, dfR], axis=1)
                
                df = pd.concat([df, df_temp], axis=0)
            except:
                print('Cannot make dataframe from {}'.format(file))


### To analyse files using metafile
else:
    
    xlfile = folder + 'metafile.xlsx'
    metafilename = folder + 'anc_metafile'
    
    metafilemaker(xlfile, metafilename, sheetname='metafile', fileformat='txt')
    
    rows, header = metafilereader(metafilename+'.txt')
    
    df=pd.DataFrame()
    
    for row in rows:
        filename=folder+row[0]
        session=int(row[1][0])
        
        sessionID = 'Rat: {}, day: {}, type: {}'.format(row[2], row[4], row[8])
        
        data = extract_data(filename, session=session, minburstlength=3)
    
    #    fig = make_histogram_fig(data[0]['lickdata'], data[1]['lickdata'], sessionID)
    #    fig.savefig('C:\\Github\\anc-for-jas\\data\\jess\\figs\\' + file + '.png')
        if row[6] == 'lr':
            df_ph1 = assemble_dataframe(data[0]['lickdata'], suffix='_ph1')
            df_ph2 = assemble_dataframe(data[1]['lickdata'], suffix='_ph2')
        else:
            df_ph1 = assemble_dataframe(data[1]['lickdata'], suffix='_ph1')
            df_ph2 = assemble_dataframe(data[0]['lickdata'], suffix='_ph2')
        
        df_id = pd.DataFrame()
        df_id['filename'] = [row[0]]
        df_id['session'] = [str(session)]
        df_id['rat'] = [row[2]]
        df_id['day'] = [row[4]]
    
        df_temp = pd.concat([df_id, df_ph1, df_ph2], axis=1)
       
        df = pd.concat([df, df_temp], axis=0)
        
### To save dataframe as excel file

writer = pd.ExcelWriter(folder + 'dataframe.xlsx')
df.to_excel(folder + 'dataframe.xlsx', 'ANC data')

pickle_out = open(folder + 'anc_dataframe.pickle', 'wb')
dill.dump(df, pickle_out)
pickle_out.close()



#        try:
#            if not self.list_of_files:
#                self.currpath = ntpath.dirname(self.filename)
#                self.list_of_files = os.listdir(self.currpath)
#            index = [x[0] for x in enumerate(self.list_of_files) if x[1] == self.shortfilename.get()]
#            newindex = index[0] + delta
#            self.filename = os.path.join(self.currpath, self.list_of_files[newindex])




