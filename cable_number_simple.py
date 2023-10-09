# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:58:37 2019

@author: Shane
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
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, LogFormatter
import seaborn as sns
import math
from scipy.spatial import distance

#import os to handle operating system
import os
#============================================================================= 
#Goal: 
#==============================================================================
#setup the data directory
datadir = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_number_data\\200416_cable_enumeration_cell_size\\"

df = pd.DataFrame()

df = pd.read_excel(datadir + '200416_cable_enumeration_cell_size.xlsx')

  
# df_1 = df[(df['d_avg'] < 6)].reset_index()  
# df_1['Bin'] = 1

# df_2 = df[(df['d_avg'] > 6) & (df['d_avg'] < 8)].reset_index()
# df_2['Bin'] = 2 

# df_3 = df[(df['d_avg'] > 8)].reset_index()
# df_3['Bin'] = 3

#=============================================================================
#parse the data into the necessary strain types for plotting
#setup df with only yBG12 cells
df_hap = pd.DataFrame()

df_hap = df.loc[(df['strain']=='yBG12')].reset_index()

#setup df with only yBG9 cells
df_dip = pd.DataFrame()

df_dip = df.loc[(df['strain']=='yBG9')].reset_index()

#setup df with only uninduced cdc28 cells
df_cdcu = pd.DataFrame()

df_cdcu = df.loc[(df['strain']=='cdc28-13_uninduced')].reset_index()

#setup df with only induced cdc28 cells
df_cdci = pd.DataFrame()

df_cdci = df.loc[(df['strain']=='cdc28-13_induced')].reset_index()


df_all = pd.DataFrame()
frames = [df_hap, df_dip, df_cdcu, df_cdci]
df_all = pd.concat(frames)
# df_all = df_all.sort_values(by='strain', ascending=False)


df_n = df.groupby([ 'n', 'strain'], sort=False).mean().reset_index()

df_n_sort = df.groupby(['n', 'strain']).mean().reset_index()

df_mean = df.groupby(['strain']).mean().reset_index()

df_std = df.groupby(['strain']).std().reset_index()

df_CI = pd.DataFrame()
df_CI['strain'] = df_mean['strain']
df_CI['mean'] = df_mean['cable_number']
df_CI['std'] = df_std['cable_number']
df_CI['interval'] = 1.96 * df_CI['std'] / np.sqrt(3) 
#=============================================================================
   
ft = 22 #font size for axis    

st = 'whitegrid' #set the style of ticks

# l = np.linspace(0, 13, 21) #lengths to plot 

ms = 200 #marker size

#set order for categorical plots
o = ['yBG12', 'yBG9','cdc28-13_uninduced','cdc28-13_induced']

#set color palette to use for plots
cmap = ["#f6511d", "#00a6ed", "#7fb800", "#ffb400"]# "#0d2c54"] 

# cmap2 = ["#f6511d", "#7fb800", "#0d2c54"] 


#=============================================================================

#plot cable number 
    
with sns.axes_style(st):
    plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
    sns.set_palette(cmap)
    
    sns.swarmplot(x='strain', y='cable_number', data = df_all, linewidth=1.0,\
                  order=o,\
                      edgecolor='k', zorder=1, size=12, dodge=True)   
        
    # ax = sns.stripplot(x='strain', y='cable_number', data = df_n, size=15,\
    #                    order=o,
    #                         color='grey', edgecolor='k', marker="s",\
    #                     linewidth=1)
    
    ax = sns.stripplot(x='strain', y='cable_number', data = df_n_sort[:4], size=15,\
                        color='grey', edgecolor='k', marker="s",\
                        linewidth=1, order = o)
        
    ax = sns.stripplot(x='strain', y='cable_number', data = df_n_sort[4:8], size=15,\
                        color='grey', edgecolor='k', marker="o",\
                        linewidth=1, order = o)
        
    ax = sns.stripplot(x='strain', y='cable_number',\
                        data = df_n_sort[8:12], size=15,\
                        color='grey', edgecolor='k', marker="^",\
                        linewidth=1, order = o) 
        
    ax = sns.pointplot(x='strain', y='cable_number', data = df_n,\
                       order=o,
                           capsize = 0.8, join=False, color='k')
        
    plt.ylabel(u'Cable number', fontsize=ft)
    plt.xlabel(None)
    ax.set(xticks=[])
    ax.tick_params('both', length=5, which='both')
    ax.yaxis.set_major_locator(ticker.MultipleLocator(4))            
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    plt.tight_layout()
    plt.ylim([0, 25])
    plt.legend([],[], frameon=False)
    # plt.savefig('230522_cable_number_categorical.svg')  
#============================================================================
    
#first parse the N means for each bin into sub-dataframes
init_hap = df_n.loc[df_n['strain']=='yBG12'] 
init_dip = df_n.loc[df_n['strain']=='yBG9'] 
init_in = df_n.loc[df_n['strain']=='cdc28-13_induced'] 
init_un = df_n.loc[df_n['strain']=='cdc28-13_uninduced'] 

#conduct t-test
stat_hap_dip, pval_hap_dip = scipy.stats.ttest_ind(init_hap['cable_number'],
                                                   init_dip['cable_number'])
#print the p-value
print(pval_hap_dip)

stat_hap_un, pval_hap_un = scipy.stats.ttest_ind(init_hap['cable_number'],
                                                   init_un['cable_number'])
#print the p-value
print(pval_hap_un)

stat_hap_in, pval_hap_in = scipy.stats.ttest_ind(init_hap['cable_number'],
                                                   init_in['cable_number'])
#print the p-value
print(pval_hap_in)

stat_dip_un, pval_dip_un = scipy.stats.ttest_ind(init_dip['cable_number'],
                                                   init_un['cable_number'])
#print the p-value
print(pval_dip_un)

stat_dip_in, pval_dip_in = scipy.stats.ttest_ind(init_dip['cable_number'],
                                                   init_in['cable_number'])
#print the p-value
print(pval_dip_in)

stat_un_in, pval_un_in = scipy.stats.ttest_ind(init_un['cable_number'],
                                                   init_in['cable_number'])
#print the p-value
print(pval_un_in)
    
    
#============================================================================
#compare quantification with published data

# with sns.axes_style(st):
#     plt.figure(figsize=(5,6))#use 3,4 for figures; 8,9 for terminal
#     sns.set_palette(cmap)
    
#     sns.swarmplot(x='strain', y='cable_number', data = df_all, linewidth=0.5,\
#                   order=['yBG12', 'cdc28-13_uninduced', 'pub_data'],\
#                       edgecolor='k', zorder=0, size=10, dodge=True)   
        
#     ax = sns.stripplot(x='strain', y='cable_number', data = df_n, size=15,\
#                        order=['yBG12', 'cdc28-13_uninduced', 'pub_data'],
#                             color='grey', edgecolor='k', marker="s",\
#                         linewidth=1)
    
        
#     ax = sns.pointplot(x='strain', y='cable_number', data = df_n,\
#                        order=['yBG12', 'cdc28-13_uninduced', 'pub_data'],
#                            capsize = 0.8, join=False, color='k')
        
#     plt.ylabel(u'Cable number', fontsize=ft)
#     plt.xlabel(None)
#     ax.set(xticks=[])
#     ax.tick_params('both', length=5, which='both')
#     ax.yaxis.set_major_locator(ticker.MultipleLocator(4))            
#     plt.rc('xtick', labelsize=ft)
#     plt.rc('ytick', labelsize=ft)
#     plt.tight_layout()
#     plt.ylim([0, 22])
#     plt.legend([],[], frameon=False) 
#     plt.savefig('210927_cable_number_categorical_pubdata.svg')    


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
        Hue to indicate replicates.        
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
    
    #calculate the standard deviation for the fit parameters
    sigma = 1.96 * (np.sqrt(np.diagonal(covar))/np.sqrt(3))
    
    #set a range of x's to plot the fit over
    x_fit = np.linspace(0,20,21)
    
    #set color palette to use for plots
    # cmap = ["#f6511d", "#00a6ed", "#7fb800", "#ffb400"] 

            
    with sns.axes_style('ticks'):
        plt.figure(figsize=(7,7))#use 3,4 for figures; 8,9 for terminal
        sns.set_palette(cmap)
    
        ax=sns.scatterplot(x, y, d, s=150, edgecolor='k', palette=c,\
                           label=None, linewidth=1.0, legend=False)
        
        ax = sns.lineplot(x_fit, powerlaw(x_fit,*pars), lw=3, color='k',\
              label= r"Fit, Slope = {0:.2f}$\pm${2:.2f}, R$^2$ = {1:.2f}"\
              .format(pars[1], r2, sigma[1]))
            
        # ax = sns.lineplot(x, powerlaw(x,1.92, 1), lw=3, color='r',\
        #       label= 'alpha=1')            
                
        ax.set(xscale="log", yscale="log")   

        for axis in [ax.xaxis, ax.yaxis]:
            formatter_min = LogFormatter(labelOnlyBase=True)
            axis.set_minor_formatter(formatter_min)
            
        plt.ylabel('Cable number', fontsize=24)    
        plt.xlabel(u'Mother cell length (${\mu}m$)', fontsize=24)
        ax.tick_params('both', length=10, which='both')
        #ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        plt.ylim([1e0, 7e1]) 
        plt.xlim([2, 2e1])
        plt.rc('xtick', labelsize=24)
        plt.rc('ytick', labelsize=24)
        plt.legend(loc='upper left', prop={'size': 18})
        plt.tight_layout()
      
        plt.savefig(sv) 
        
    
#==============================================================================
#make plots for bud/mom, total_cell_corr_int, cable_corr_int, patch_corr_int
#colors to use: '#f6511d', '#00a6ed', '#7fb800', '#ffb400', '#1CCAD8'
             
scaling_plot(powerlaw, cof_det,\
             df_all['d1'], df_all['cable_number'], df_all['strain'],\
             cmap,\
                 '230522_cable_number_scaling.svg')
    
#==============================================================================
#make a plot of total intensity vs volume with the bud cell and mother cells
#in different colors
#setup df with only mother cdc28 cells

# with sns.axes_style('ticks'):
#     f, ax = plt.subplots(figsize=(8, 8))
#     # plt.loglog(df['d_avg'],powerlaw(df['d_avg'],*pars_v),'k--', \
#     #             label="Fit, Slope = " + str(pars_v[1].round(2)))

#     plt.loglog(df_bud['volume'],df_bud['total_cell_corr_int'], \
#                 marker='o', markersize=13, linestyle='None',
#                 mew=1, mfc='#00a6ed', mec='k', label='Bud, cell')
        
#     plt.loglog(df_mom['volume'],df_mom['total_cell_corr_int'], \
#                 marker='o', markersize=13, linestyle='None',
#                 mew=1, mfc='#f6511d', mec='k', label='Mom, cell')
#     # ax = sns.scatterplot(x=df['volume'], y=df['total_cell_corr_int'],\
#     #                 edgecolor='k')
        
        
    
#     plt.ylabel('log(Integrated density)', fontsize=24)    
#     plt.xlabel('log(Cell volume)', fontsize=24)
#     plt.legend(loc='upper left')
#     plt.rc('xtick', labelsize=24)
#     plt.rc('ytick', labelsize=24)
#     plt.tight_layout()
#     plt.savefig('210325_Bud_Mom_totalcellint_diff_colors.svg')    

