# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 14:44:13 2021

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
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import seaborn as sns
import math
from math import sqrt
from scipy.spatial import distance
from scipy.signal import chirp, find_peaks, peak_widths
from sklearn import metrics


#import os to handle operating system
import os
#=============================================================================
#Goal: Import phalloidin stained actin cable traces analyzed using 
#'220621_actin_intensity_profile_calcs.py' and conduct fitting and plotting 
#of this data.
#=============================================================================
#Define the directory containing the anaylzed and compiled data.
datadir = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_intensity_profile_data\\final_data\\"

#=============================================================================
#Initialize an empty dataframe to write the data to.
df = pd.DataFrame()
#Read in the data to the empty dataframe.
df = pd.read_csv(datadir + '230510_cdc28_all_cable_traces_rescaled.csv')
#remove all of the data in bin 2
df = df[df.bin != 2]
#output cleaned data
df.to_csv(datadir + '230510_cdc28_all_cable_traces_rescaled_clean.csv', index=False)
#Determine the length of each cable analyzed (i.e. the maximum length measured)
df['cable_length'] = df.groupby('cell')['Distance_(microns)'].transform('max')
#initialize dataframe to write cell-specific data to
df_cell = pd.DataFrame()
df_cell = df.groupby(['cell']).max().reset_index()
#Calculate the mean for each experimental replicate in each bin.
df_all_mean = df.groupby(['bin', 'n'], sort=True).mean().reset_index()

#Calculate the 95% confidence intervals for each bin
df_mean = df.groupby(['bin']).mean().reset_index()

df_std = df.groupby(['bin']).std().reset_index()

df_CI = pd.DataFrame()
df_CI['bin'] = df_mean['bin']
df_CI['mean'] = df_mean['cable_length']
df_CI['std'] = df_std['cable_length']
df_CI['interval'] = 1.96 * df_CI['std'] / np.sqrt(3)

#=============================================================================
#Define the directory containing the anaylzed and compiled data.
datadir3 = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_intensity_profile_data\\final_data\\"

#=============================================================================
#Initialize an empty dataframe to write the data to.
df_wt = pd.DataFrame()
#Read in the data to the empty dataframe.
df_wt = pd.read_csv(datadir3 + '230421_all_cable_traces_rescaled_trim.csv')
#Determine the length of each cable analyzed (i.e. the maximum length measured)
df_wt['cable_length'] = df_wt.groupby('cell')['Distance_(microns)'].transform('max')
#initialize dataframe to write cell-specific data to
df_cell_wt = pd.DataFrame()
df_cell_wt = df_wt.groupby(['cell']).max().reset_index()
#Calculate the mean for each experimental replicate in each bin.
df_all_mean_wt = df_wt.groupby(['bin', 'n'], sort=True).mean().reset_index()
#normalize the relative intenstiy values between 0-1

#=============================================================================
#Initalize the plotting parameters.   
ft = 24 #font size for axis    
# ft2 = 18 #font size for axis    

st = 'whitegrid' #set the style of ticks
# l = np.linspace(0, 13, 21) #lengths to plot 
ms = 200 #marker size
#set color palette to use for plots
cmap = ["#7fb800", "#ffb400"]#, "#0d2c54","#f6511d", "#00a6ed"]
#=============================================================================
#plot the cable lengths within each bin

with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    
    sns.set_palette(cmap)  
    
    ax = sns.swarmplot(x= df_cell['bin'], y= df_cell['cable_length'], \
                        linewidth=0.7, edgecolor='k', zorder=1, \
                            size=12, dodge=True)

    ax = sns.stripplot(x='bin', y='cable_length', hue='n', \
                        data = df_all_mean, size=15,\
                        color='grey', edgecolor='k', marker="s",\
                        linewidth=1) 
        
    ax = sns.pointplot(x='bin', y='cable_length', data = df_all_mean,\
                        capsize = 0.8, join=False, color='k')
        

    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.xlabel(u'Bin', fontsize=ft)
    
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2)) 
     
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    # plt.xlim([0, 13])
    plt.ylim([0, 13])
    plt.legend([],[], frameon=False)
    plt.tight_layout()
    # plt.savefig('230510_actin_cable_int_profile_cable_len.svg') 

#=============================================================================
#plot the cell lengths within each bin

with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    
    sns.set_palette(cmap) 
    
    ax = sns.swarmplot(x= df_cell['bin'], y= df_cell['cell_length'], \
                        linewidth=0.5, edgecolor='k', zorder=1, \
                            size=12, dodge=True)

    ax = sns.stripplot(x='bin', y='cell_length', hue='n', \
                        data = df_all_mean, size=15,\
                        color='grey', edgecolor='k', marker="s",\
                        linewidth=1) 
        
    ax = sns.pointplot(x='bin', y='cell_length', data = df_all_mean,\
                        capsize = 0.8, join=False, color='k')
        

    plt.ylabel(u'Cell length (${\mu}m$)', fontsize=ft)
    plt.xlabel(u'Bin', fontsize=ft)
    
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2)) 
      
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    # plt.xlim([0, 13])
    plt.ylim([0, 16])
    
    plt.legend([],[], frameon=False)
    plt.tight_layout()
    # plt.savefig('230510_actin_cable_int_profile_cell_len.svg')  

#=============================================================================    
#plot the initial intensity of cables closest to the bud neck
  
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    
    sns.set_palette(cmap)  
    
    ax = sns.swarmplot(x= df_cell['bin'], y= df_cell['int_0'], linewidth=0.7,\
                  edgecolor='k', zorder=1, size=12, dodge=True, label=None)

    ax = sns.stripplot(x='bin', y='int_0', hue='n', data = df_all_mean, size=15,\
                        color='grey', edgecolor='k', marker="s",\
                        linewidth=1, label=None) 
        
    ax = sns.pointplot(x='bin', y='int_0', data = df_all_mean,\
                        capsize = 0.8, join=False, color='k', label='Bin')
        

    plt.ylabel(u'Cable intensity at bud neck (AU)', fontsize=ft)
    plt.xlabel(u'Bin', fontsize=ft)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.major.formatter._useMathText = True     
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    # plt.xlim([0, 13])
    plt.ylim([0, 8e4])
    
    plt.legend([],[], frameon=False)
    # plt.savefig('230510_actin_cable_int_profile_initial_int.svg') 

#==============================================================================    
#Run a t-test on the inital intensity values

#first parse the N means for each bin into sub-dataframes
intit_bin1 = df_all_mean.loc[df_all_mean['bin']==1.0] 
intit_bin3 = df_all_mean.loc[df_all_mean['bin']==3.0] 
#conduct t-test for initial intensity
stat_int, pval_int = scipy.stats.ttest_ind(intit_bin1['int_0'],
                                           intit_bin3['int_0'])
#print the p-value
print(pval_int)

#conduct t-test for cable length
stat_len, pval_len = scipy.stats.ttest_ind(intit_bin1['cable_length'],
                                           intit_bin3['cable_length'])
#print the p-value
print(pval_len)

#conduct t-test for cell length
stat_celllen, pval_celllen = scipy.stats.ttest_ind(intit_bin1['cell_length'],
                                           intit_bin3['cell_length'])
#print the p-value
print(pval_celllen)

#==============================================================================    
#import curve_fit from scipy
from scipy.optimize import curve_fit

#write a function to fit using a single exponential with curve fit
def exp_fit (t, a, k, c):
    '''
    Parameters
    ----------
    t : independent varibale (i.e. time or distance).
    a : peak amplitude (i.e. y-intercept).
    k : rate constant (either decay or time constant).
    c : offset (constant value)

    Returns
    -------
    y_fit : function fit to the data
    
    '''
    y_fit = a * np.exp(-t / k) + c
    return y_fit


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
    
    res = y - exp_fit(x, *z)
    
    ss_res = np.sum(res**2)
    
    ss_tot = np.sum(( y - np.mean(y))**2)

                    
    r_squared = 1 - (ss_res / ss_tot)

    return (r_squared)

#==============================================================================
#write a function to compute the fold change in values and 
#propogate the errors for this comparison

def fold_change_err(mean_1, std_1, mean_2, std_2):
    '''
    Parameters
    ----------
    mean_1: the mean of quantity one
    std_1: the standard deviation of quantity one
    mean_2: the mean of quantity two
    std_2: the standard deviation of quantity two

    Returns
    -------
    fc : fold change or ratio of mean_1:mean_2
    err: the error associate with the fold_change
    '''
    fc = mean_1 / mean_2
    
    err = fc * np.sqrt((std_1/mean_1)**2 + (std_2/mean_2)**2)
    
    return fc, err
    
cell_length_fc = fold_change_err(intit_bin3['cell_length'].mean(),
                                 intit_bin3['cell_length'].std(),
                                 intit_bin1['cell_length'].mean(),
                                 intit_bin1['cell_length'].std())

cable_length_fc = fold_change_err(intit_bin3['cable_length'].mean(),
                                 intit_bin3['cable_length'].std(),
                                 intit_bin1['cable_length'].mean(),
                                 intit_bin1['cable_length'].std())
    
#==============================================================================
#Fit the data in bin 1 to a single exponential and calculate the coefficient
#of determination

#group all of the data from bin 1 into a single dataframe
bin1 = df[(df['bin'] == 1)].reset_index()  
bin1 = bin1.groupby(['n', 'Distance_(microns)']).mean().reset_index()
#calculate the number of cells in bin 1
# n1  = bin1.groupby(['cell']).count()

#normalize the relative intenstiy values between 0-1
import sklearn
from sklearn import preprocessing

#fit the cell length normalized data from bin1 to a single exponential 
pars_norm1, covar_norm1 =  curve_fit(exp_fit, bin1['cell_len_norm'], \
                                     bin1['int_norm_rescale'], \
                              absolute_sigma = False, p0=(1, 0.3, 1))
#calculate the coefficient of determination for the fit    
r2_norm1 = cof_det(bin1['int_norm_rescale'], bin1['cell_len_norm'], pars_norm1)
#calculate the standard deviation for the fit
sigma_norm1 = 1.96 * (np.sqrt(np.diagonal(covar_norm1))/np.sqrt(5))  

#do the same for the raw data
pars_dist1, covar_dist1 =  curve_fit(exp_fit, bin1['Distance_(microns)'],\
                              bin1['int_norm_rescale'], \
                              absolute_sigma = False, p0=(1, 1, 1))  
#calculate the coefficient of determination for the fit        
r2_dist1 = cof_det(bin1['int_norm_rescale'], bin1['Distance_(microns)'], pars_dist1)
#calculate the standard deviation for the fit
dist_1_std = np.sqrt(np.diagonal(covar_dist1))
sigma_dist1 = 1.96 * (dist_1_std/np.sqrt(5)) 

#==============================================================================
#Fit the data in bin 3 to a single exponential and calculate the coefficient
#of determination

#group all of the data from bin 3 into a single dataframe
bin3 = df[(df['bin'] == 3)].reset_index()  
bin3 = bin3.groupby(['n', 'Distance_(microns)']).mean().reset_index()

#fit the cell length normalized data from bin3 to a single exponential 
pars_norm3, covar_norm3 =  curve_fit(exp_fit, bin3['cell_len_norm'],\
                              bin3['int_norm_rescale'], \
                              absolute_sigma = False, p0=(1, 0.3, 1))  
    
#calculate the coefficient of determination for the fit        
r2_norm3 = cof_det(bin3['int_norm_rescale'], bin3['cell_len_norm'], pars_norm3)
#calculate the standard deviation for the fit
sigma_norm3 = 1.96 * (np.sqrt(np.diagonal(covar_norm3))/np.sqrt(5)) 

#do the same for the raw data
pars_dist3, covar_dist3 =  curve_fit(exp_fit, bin3['Distance_(microns)'],\
                              bin3['int_norm_rescale'], \
                              absolute_sigma = False, p0=(1, 1, 1))  
    
#calculate the coefficient of determination for the fit        
r2_dist3 = cof_det(bin3['int_norm_rescale'], bin3['Distance_(microns)'], pars_dist3)
#calculate the standard deviation for the fit
dist_3_std = np.sqrt(np.diagonal(covar_dist3))
sigma_dist3 = 1.96 * (dist_3_std/np.sqrt(5)) 
#=============================================================================
#compute the fold_change in lambda
lambda_fc = fold_change_err(pars_dist3[1], dist_3_std[1],
                            pars_dist1[1], dist_1_std[1])
print(lambda_fc)

#=============================================================================
#group all of the data from bin 1 into a single dataframe
wt_mean = df_wt.groupby(['n', 'Distance_(microns)']).mean().reset_index()

#do the same for the wt data
pars_dist_wt, covar_dist_wt =  curve_fit(exp_fit, wt_mean['Distance_(microns)'].round(1),\
                              wt_mean['int_norm_rescale'], \
                              absolute_sigma = False, p0=(1, 1,1))  
#calculate the coefficient of determination for the fit        
r2_dist_wt = cof_det(wt_mean['int_norm_rescale'], wt_mean['Distance_(microns)'].round(1), pars_dist_wt)
#calculate the standard deviation for the fit
sigma_dist_wt = 1.96 * (np.sqrt(np.diagonal(covar_dist_wt))/np.sqrt(3)) 
 
#==============================================================================
#define the values to plot the fitted parameters over   
l_fit_d = np.linspace(0, 10, 101)
l_fit_n = np.linspace(0, 1, 101)
#==============================================================================
#Plot the normalized intensity versus the normalized cell lengths
with sns.axes_style(st):
    plt.figure(figsize=(6,6))#use 3,4 for figures; 8,9 for terminal
    
    sns.color_palette("flare", as_cmap=True)    
     
    ax = sns.lineplot(x= bin1['cell_len_norm'].round(2),
                      y= bin1['int_norm_rescale'],\
                          lw=4, label=None, color='#7fb800')
            
    ax = sns.lineplot(x= bin3['cell_len_norm'].round(2),
                      y= bin3['int_norm_rescale'],\
                          lw=4, label=None, color='#ffb400') 
        
    sns.lineplot(l_fit_n, exp_fit(l_fit_n, *pars_norm1), color = 'k', \
                  label= r"Decay length = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                    .format(pars_norm1[1],r2_norm1,sigma_norm1[1]), lw=2,
                    linestyle='--')         

    sns.lineplot(l_fit_n, exp_fit(l_fit_n, *pars_norm3), color = 'grey', \
                  label= r"Decay length = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                    .format(pars_norm3[1],r2_norm3,sigma_norm3[1]), lw=2,
                    linestyle='-')         
                         
    plt.ylabel(u'Relative Cable Intensity', fontsize=ft)
    plt.xlabel(u'Cable length / Mother cell length', fontsize=ft)
    
    ax.tick_params('both', length=5, which='both')
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(2)) 
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(2)) 
          
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.xlim([0, 1])
    plt.ylim([0, 1.1])
    
    plt.legend(loc='upper right', prop={'size': 14})
    # plt.savefig('230510_actin_cable_int_profile_cdc28_norm.svg')  
   
#==============================================================================
#Plot the normalized intensity vs cable length 
    
with sns.axes_style(st):
    plt.figure(figsize=(6,6))#use 3,4 for figures; 8,9 for terminal
    
    sns.color_palette("flare", as_cmap=True)    
            
    ax = sns.lineplot(x= bin1['Distance_(microns)'].round(1),
                      y= bin1['int_norm_rescale'], lw=4, color='#7fb800')
    
    ax = sns.lineplot(x= bin3['Distance_(microns)'].round(1), 
                      y= bin3['int_norm_rescale'], lw=4, color='#ffb400')        
        
    sns.lineplot(l_fit_d, exp_fit(l_fit_d, *pars_dist1), color = 'k', \
                  label= r"Decay length = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                    .format(pars_dist1[1],r2_dist1,sigma_dist1[1]), lw=2,
                    linestyle='--') 

    sns.lineplot(l_fit_d, exp_fit(l_fit_d, *pars_dist3), color = 'grey', \
                  label= r"Decay length = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                    .format(pars_dist3[1],r2_dist3,sigma_dist3[1]), lw=2,
                    linestyle='-') 
    
    plt.ylabel(u'Relative Cable Intensity', fontsize=ft)
    plt.xlabel(u'Cable length (${\mu}m$)', fontsize=ft)
    
    ax.tick_params('both', length=5, which='both')
          
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.xlim([0, 10])
    plt.ylim([0, 1.1])
    
    plt.legend(loc='upper right', prop={'size': 14})
    # plt.savefig('230510_actin_cable_int_profile_cdc28_dist.svg')
    
#==============================================================================
#Plot the normalized intensity versus the measured cell lengths
    
    
with sns.axes_style(st):
    plt.figure(figsize=(6,6))#use 3,4 for figures; 8,9 for terminal         

    ax = sns.lineplot(x= df_wt['Distance_(microns)'].round(1),
                      y= df_wt['int_norm_rescale'],lw=4, color='#C200FB',
                      label=None)
                   
    sns.lineplot(l_fit_d, exp_fit(l_fit_d, *pars_dist_wt),
                 color = 'k', \
                  label= r"Decay length = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
                    .format(pars_dist_wt[1],r2_dist_wt,sigma_dist_wt[1]), lw=2,
                    linestyle='--')
        
                         
    plt.ylabel(u'Relative Cable Intensity', fontsize=ft)
    plt.xlabel(u'Cable length (${\mu}m$)', fontsize=ft)
    
    ax.tick_params('both', length=5, which='both')
          
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.xlim([0, 6])
    plt.ylim([0, 1.1])
    
    plt.legend(loc='upper right', prop={'size': 14})
    # plt.savefig('220510_actin_cable_int_profile_wt.svg')      
    
 
    
    

    
