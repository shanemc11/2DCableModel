# -*- coding: utf-8 -*-
"""
Created on Fri May 14 12:51:29 2021

@author: shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import ScalarFormatter, LogFormatter
import seaborn as sns
import math
from scipy.spatial import distance

#import os to handle operating system
import os
#=============================================================================
#Goal: Compare measured cable architecture parameters from different strains/
#mutants. 

#This script will import the data file from the indicated data directory 
#and fit these data using the power law, determine the coefficient of 
#determination from this fit, and plot the raw data onto a loglog plot with
#the results of the above analyses. The resulting plot will be saved as 
#a svg file in the current directory.
#=============================================================================
#indicate the data directory
datadir = "insert data directory here"

#initialize empty dataframes
df = pd.DataFrame()
   
df = pd.read_csv(datadir + 'filename') 
#=============================================================================
#import curve_fit from scipy
from scipy.optimize import curve_fit

#write a function to fit using the power law with curve fit
def powerlaw(x,a,b):
    '''
    Parameters
    ----------
    x : The size dimension to scale intensity/length data with.
    a : Constant.
    b : Scaling exponent.

    Returns
    -------
    Use to fit data and extract scaling exponent for various cellular
    dimensions.

    '''
    y = a*(x**b)
    
    return y

#==============================================================================
#write a function to calculate the coefficient of determination for powerlaw
#fits

def cof_det(y, x, z):
    '''
    Parameters
    ----------
    y : dependent variable from powerlaw fits (ex: intensity, length).
    x : independent varibale from powerlaw fits (cells size dimensions).
    z : fitting parameter from powerlaw fits (scaling exponent and constant).

    Returns
    -------
    r_squared : coefficient of determination of powerlaw fit.

    '''
    
    res = y - powerlaw(x, *z)
    
    ss_res = np.sum(res**2)
    
    ss_tot = np.sum(( y - np.mean(y))**2)

                    
    r_squared = 1 - (ss_res / ss_tot)

    return (r_squared)

#==============================================================================
#write a function to calculate the powerlaw scaling parameters, the 
#coefficient of determination, and plot the results on a loglog plot

def scaling_plot(powerlaw, cof_det, x, y, d, c, sv):
    '''

    Parameters
    ----------
    powerlaw : function
        Fits the data with the powerlaw to measure the scaling coeffcient and
        constant.
    cof_det : function
        Calculates the coefficient of determination.
    x : variable
        Dependent variable to use in the above functions.
    y : variable
        Independent variable to use in the above functions.
    d : variable
        Cell length measured from the bud neck to cell rear.        
    c : string
        Color for markers.        
    txt : string
        Label for plot.
    sv : string
        Label for plot file during save.        
                

    Returns
    -------
    Results from fitting the data with the powerlaw and a loglog plot.

    '''
    pars, covar = curve_fit(powerlaw,x,y)
    
    r2 = cof_det(y, x, pars)
    
    #set the range of lengths to plot the fit
    l = np.linspace(3, 13, 20)
    
    #calculate the standard deviation for the fit parameters
    sigma = 1.96 * (np.sqrt(np.diagonal(covar)))/np.sqrt(3)
            
    with sns.axes_style('ticks'):
        plt.figure(figsize=(7,7))#use 3,4 for figures; 8,9 for terminal
        # sns.set_palette(cmap)
    
        ax=sns.scatterplot(x=x, y=y, s=100, edgecolor='k', color=c,\
                           label=None, linewidth=1.0)
        
        ax = sns.lineplot(x=l, y=powerlaw(l,*pars), lw=3, color='k',\
              label= r"Fit, Slope = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
              .format(pars[1], r2, sigma[1]))
            
        # ax = sns.lineplot(x, powerlaw(x,0.006, 1), lw=3, color='r',\
        #       label= 'alpha=2')            
                
        ax.set(xscale="log", yscale="log")   

        for axis in [ax.xaxis, ax.yaxis]:
            formatter_min = LogFormatter(labelOnlyBase=True)
            axis.set_minor_formatter(formatter_min)
            
        plt.ylabel('Bnr1-GFP$^{Envy}$ Intensity', fontsize=24)    
        plt.xlabel(u'Mother cell length (${\mu}m$)', fontsize=24)
        ax.tick_params('both', length=10, which='both')
        #ax.yaxis.set_major_locator(ticker.MultipleLocator(2))            
        plt.rc('xtick', labelsize=24)
        plt.rc('ytick', labelsize=24)
        plt.xlim([3, 13])
        plt.ylim([1e4, 2e6])
        plt.legend(loc='upper left', prop={'size': 20})
        plt.tight_layout()
      
        plt.savefig(sv) 

#==============================================================================
#colors to use: '#f6511d', '#00a6ed', '#7fb800', '#ffb400', '#1CCAD8'
             
# scaling_plot(powerlaw, cof_det,\
#               df['cell_diameter'],df['int_corr'],df['n'], \
#               '#F28D35', '#outputname.svg')  
    
   