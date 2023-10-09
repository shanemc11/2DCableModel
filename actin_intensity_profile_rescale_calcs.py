# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:21:47 2022

@author: shane
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 19:35:04 2021

@author: shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
import operator
from operator import truediv
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import math
from math import sqrt
from scipy.spatial import distance
from scipy.signal import chirp, find_peaks, peak_widths
from sklearn import metrics
import sklearn
from sklearn import preprocessing


#import os to handle operating system
import os
#=============================================================================
#Goal: Import flourescence profiles of phalloidin stained actin cables and 
# calculate background subtracted intensities. Also, correlate these intensity
# values with cell size parameters and sort cables into bins based on these
# cell size parameters. This script will read in intensity profiles measured
# in FIJI/ImageJ that are saved as csv files and correlate these measurements
# with a excel files containing background intensity and cell size parameters.
#=============================================================================
#Indicate the directory containing the folders with cable traces and the
# file containing the cell dimension data.

datadir = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_intensity_profile_data\\"

#==============================================================================
#Write a function to compile these traces into a single dataframe to perform
# calcualtions and normalization of the data.

def profile_calcs(datadir, trace_dir, cell_dims, s, output1):
    '''
    

    Parameters
    ----------
    datadir : directory
        Directory to folder with files for analysis.
    trace_dir: sub directory
        Indicates the folder name that contains the cable intensity traces.
    cell_dims: str
        Name of csv file that has list of cell name (same as profile file 
        name) and length of the cell in a column named 'length'
    s: int
        Size cutoff for bin assignment.        
    output : str
        Output file name.

    Returns
    -------
    Calculations for analysis of intensity profiles of F-actin in cells
    of different size.

    '''
    #initialize empty dataframes to write the data to
    df = pd.DataFrame()
    df_calcs = pd.DataFrame()
    celldims = pd.DataFrame()
    
    


    for f in glob.glob(datadir + trace_dir + '*' + '.csv'):
    
        #open each file
        df = pd.read_csv(f)
        
        #read the cell dimension data into a new dataframe
        celldims = pd.read_excel(datadir + cell_dims)
            
        #set the initial intensity value as the first value for each trace
        df['int_0'] = df['Gray_Value'][0]
        
        #find the max intensity value from each trace
        df['int_max'] = df['Gray_Value'].max()
        
        #set the total length of the cable as the length values measured in
        #the trace file - this simplifies calling length parameters downstream                
        df['cable'] = df['Distance_(microns)']
       
        #assign the file name as the identifier for the cell that was measured
        df['cell'] = os.path.basename(f)
        
        #write the background intensity for each image to the corresponding cell
        df['bckgrd'] =  df.cell.map(celldims.set_index('cell').bckgrd)
        
            #subtract the background intensity from each cell       
        df['int_corr'] = df['Gray_Value'] - df['bckgrd']
        
        #normalize the measure intensity values to the initial intensity
        df['int_norm'] = df['int_corr'] / (df['int_max'] - df['bckgrd'])
        
        df['int_norm_rescale'] = sklearn.preprocessing.minmax_scale(df['int_norm'], feature_range=(0, 1), axis=0, copy=True)
        
            #normalize the cable length by the cell length for each cell
        df['cell_len_norm'] =  df.cable.div(df.cell.map(celldims.set_index('cell').diameter))
        
        #write the cell length from the cell_dims file to the dataframe
        df['cell_length'] =  df.cell.map(celldims.set_index('cell').diameter)
    
        
        #write the experimental N to the file for downstream stats    
        df['n'] = df.cell.map(celldims.set_index('cell').n)
        
        #assign cells to a bin based on cell length
        df['bin'] = df['cell_length']
        df.loc[df['cell_length'] <= s, 'bin'] = 1
        df.loc[(df['cell_length'] > s) & 
                     (df['cell_length'] < 7.9), 'bin'] = 2
        df.loc[df['cell_length'] > 7.9, 'bin'] = 3

        #append data from all traces into a single dataframe
        df_calcs = df_calcs.append(df)
    

    df_calcs.to_csv(datadir + output1, index=False)
        
    return df_calcs
    
    
data = profile_calcs(datadir, 'all_cable_traces\\', '210916_cell_dimensions.xlsx',\
                     6.1, '230510_cdc28_all_cable_traces_rescaled.csv')

    
