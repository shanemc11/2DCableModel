# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:23:32 2020

@author: Shane
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
from scipy.signal import find_peaks

#import os to handle operating system
import os
#=============================================================================
#Goal: Import flourescence profiles of phalloidin stained cells and detect 
# peaks of flourescence that correspond to cables. Enuermate peaks for 
# each cell and wrtie to a csv file for plotting.
#=============================================================================

# datadir = "#indicate directory containing files to analyze"

#=============================================================================
#initalize data frame to append all data 
df = pd.DataFrame()

df4 = pd.DataFrame()

t = 0.2 #peak height threshold; set to 1 for no threshold

#Read in each trace from the directory
for f in glob.glob(datadir + '*' + '.csv'):
    
    #open each file
    df = pd.read_csv(f)
    
    #set the height threshold for each trace; 
    h = ((df['Y'].max() - df['Y'].min()) * t) + df['Y'].min()
    
    #trace peaks for each trace
    peaks, _ = find_peaks(df['Y'], height = h)
    
    #initialize an empty dataframe to write detected peaks to
    df2 = pd.DataFrame()
    
    #concatenate the x-coordinates for each cell
    df2 = pd.concat([df2, pd.DataFrame(peaks)], axis=1)
    
    #initialize a dataframe to write the detected peaks to
    df3 = pd.DataFrame()
    
    #write the detected cables to the new dataframe and name the column
    df3['cable_number'] = df2[0]
    
    #indicate where the file where cables were detected
    df3['cell'] = os.path.basename(f)
    
    #append data from all traces into a single dataframe
    df4 = df4.append(df3)

#initialize a final dataframe for counting the detected peaks
df5 = pd.DataFrame()
#count the number of peaks for each cell
df5 = df4.groupby(['cell']).count().reset_index()


# #output the data to csv
# outputdir = '\#indicate directory for data output\'
# df5.to_csv(datadir + '#filename',\
#             index=False)



