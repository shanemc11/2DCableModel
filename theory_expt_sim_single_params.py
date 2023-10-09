# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 12:45:44 2023

@author: shane
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker

import random as rand
from statsmodels.distributions.empirical_distribution import ECDF
import seaborn as sns 
sns.set()
#import curve_fit from scipy
from scipy.optimize import curve_fit

#=============================================================================
#Goal: This version of the simulation simulates the addition of filaments
#whose lengths are exponential distributed. This is different from prior
#version which used a fixed value for filament length.

#=============================================================================
#These functions are used to simulate cable trajectories and 
#length distributions

def CableEvol(prob, TTime, f_len):
    '''
    This function simulates the assembly of a single cable lane given the 
    rules of the treadmilling model.
    
    Parameters
    ----------
    prob : float
        Probability of filament in the cable bundle being removed.
    TTime : int
        Time to run the simulation (seconds).
    f_len : float
        Length of filaments (subunits) that are added and removed from the
        cable bundle.
        
    Returns
    -------
    LCable : Array
        Returns an array containing the length (number of filaments) of a 
        single cable 'lane' as a function of simulation time.

    '''
    p = prob #probability filament is removed
    TotalTime = TTime #time to run simulation
    LCable = np.zeros(TotalTime) #array of cable lengths over time (output)
    Cable = np.zeros(TotalTime) #current configuration of the cable.
                                #1 in position i indicates the presence
                                #of a filament; 0 is the absence.
    Cable[0] = 0 #initial cable length 
    LCable[0] = 0 #initial cable length

    #iterate through the time range specified
    for time in range(TotalTime):
        #at each time step insert a single filament
        Cable = np.insert(Cable,0,1)
        #measure the number of filaments in the cable lane
        Len = Cable.size
        #iterate through each position in the cable lane
        for pos in range(1,Len):
            #if the position is occupied by a filament roll the dice to
            #determine if that filament remains (stays 1) or is removed 
            # (converted to zero)             
            if Cable[pos] == 1:
                dice = rand.random() 
                if dice < p:
                    Cable[pos] = 0 
        #trim any trailing zeros from the lane
        Cable = np.trim_zeros(Cable)
        #measure the number of filaments in the lane at each time step
        LCable[time] = Cable.size
    return LCable


def CompCableEvol(p, TTime, NoOfCables, f_len):
    '''
    This function simulates the assembly of a multiples cable lanes given the 
    rules of the treadmilling model.

    Parameters
    ----------
    prob : float
        Probability of filament in the cable bundle being removed.
    TTime : int
        Time to run the simulation (seconds).
    NoOfCables : int
        Number of formins that are assembling independent cable lanes.        
    f_len : float
        Length of filaments (subunits) that are added and removed from the
        cable bundle.

    Returns
    -------
    LMax : Array
        Returns an array containing the length (number of filaments) of the 
        longest cable 'lane' as a function of simulation time.

    '''
    prob = p #probability filament is removed
    Time = TTime #time to run simulation
    NoCables=NoOfCables #number of cable lanes in the bundle
    #run the simulation for a single cable lane
    Lcables = np.array([CableEvol(prob,Time, f_len)])
    #iterate of the number of lanes in the bundle
    for i in range (1,NoCables):
        Lcables = np.concatenate((Lcables,[CableEvol(prob,Time, f_len)]))
    #initialize an empty array    
    LMax=Lcables[0]
    #iterate through the cable lanes and record the max length at each time
    #step of the simulation
    for i in range (1, NoCables):
        LMax = np.maximum(LMax, Lcables[i])
 
    return LMax

#=============================================================================
#These functions are used to simulate cable tapering.

def FilaEvol(prob, TTime):
    '''
    This function evolves the cable length according to the rules for a 
    time TTime and returns a list of Cable lengths. This is similar to the 
    CableEvol function, but is adapted to permit the tracking of filament 
    density across the lane for comparisons with tapering cables.

    Parameters
    ----------
    prob : float
        Probability of filament in the cable bundle being removed.
    TTime : int
        Time to run the simulation (seconds).

    Returns
    -------
    Cable : Array
        Returns an array containing the occupancy of each position in a 
        single cable 'lane' as a function of simulation time.


    '''
    p = prob #probability filament is removed
    TotalTime = TTime #time to run simulation
    LCable = np.zeros(TotalTime) #array of cable lengths over time (output)
    Cable = np.zeros(TotalTime) #current configuration of the cable.
                                #1 in position i indicates the presence
                                #of a filament; 0 is the absence.
    Cable[0] = 0 #initial cable length 
    LCable[0] = 0 #initial cable length
    #iterate through the time range specified    
    for time in range(TotalTime):
        #at each time step insert a single filament       
        Cable = np.insert(Cable,0,1)
        #measure the number of filaments in the cable lane        
        Len = Cable.size
        #iterate through each position in the cable lane        
        for pos in range(1,Len):
            #if the position is occupied by a filament roll the dice to
            #determine if that filament remains (stays 1) or is removed 
            # (converted to zero)             
            if Cable[pos] == 1:
                dice = rand.random() 
                if dice < p:
                    Cable[pos] = 0
        #trim any trailing zeros from the lane                    
        Cable = np.trim_zeros(Cable)
        #record which positions are occupied and which are empty in the
        #cable lane - This is the step the makes this function different
        #from the CableEvol function
        Cable = Cable
    return Cable

def FilaCompEvol(p, TTime, NoOfCables):
    
    '''
    This function iterates the FilaEvol function for the number of specified
    formins

    Parameters
    ----------
    prob : float
        Probability of filament in the cable bundle being removed.
    TTime : int
        Time to run the simulation (seconds).
    NoOfCables : int
        Number of formins that are assembling independent cable lanes.        

    Returns
    -------
    all_filaments : DataFrame
        Returns an array containing the occupancy of each position in a 
        single cable 'lane' as a function of simulation time.

    '''
    prob = p #probability filament is removed
    Time = TTime #time to run simulation
    NoCables=NoOfCables #number of cable lanes in the bundle
  
    Lfilaments = pd.DataFrame() #initialize an empty dataframe
    all_filaments = pd.DataFrame() #initialize an empty dataframe
    #iterate of the number of lanes in the bundle      
    for i in range (1,NoCables):
        #record the final length of each lane and the occupancy of that lane
        Lfilaments['length'] = ([FilaEvol(prob,Time)])
        #indicate the lane number
        Lfilaments['cable'] = i
        #append the results for all lane together
        all_filaments = all_filaments.append(Lfilaments)

    return all_filaments
#=============================================================================
#These functions are used to compare a single parameter with the analytic
#results and experimental measurements.

def cable_sim(kval, time, f_len, num_f, k_plus, TotalTraj):
    '''
    Parameters
    ----------
    kval : Float
        Range of k-minus rates (per second) to use in simulation.
    time : Int
        Number of time steps to run the simulation.
    f_len : Int
        Filament length (subunits) to use in simulation.
    num_f : Int
        Formin number to use in simulation.
    k_plus : Float
        Rate of filament assembly (subunits/sec).
    TotalTraj : Int
        Number of independent simulations to run.

    Returns
    -------
    data_accum : DataFrame
        Returns an dataframe containing the maximum cable length at each time
        step and the parameters used during that simulation.

    '''
    #initalize an empty dataframe to store the results of the simulation
    data_accum = pd.DataFrame()
    #iterate over the number of independent simulations to run
    for i in range(TotalTraj):
        #compute the probability of filament removal for each
        #parameter set
        prob_var =   1 - np.exp(-1*(kval/(k_plus/f_len)))
        #initialize an empty dataframe
        data = pd.DataFrame()
        #run simulation        
        Lcable = CompCableEvol(prob_var, time, num_f, f_len)
        #convert array values to lengths (microns)        
        data['length'] = Lcable * (f_len*2.7/1000)
        #rescale lengths to start growth at zero length
        data['length_rescale'] = (data['length']-data['length'][0])
        #scale the time into seconds
        data['time'] = np.arange(0, time, 1)*(f_len)/k_plus
        #record the parameters used in each simulation run        
        data['k_minus'] = kval * f_len
        data['iteration'] = i+1
        #append the data for each iteration        
        data_accum = data_accum.append(data, ignore_index=True)
    return data_accum

def TaperSim (kval, time, f_len, num_f, k_plus , TotalCables):
    '''
    Parameters
    ----------
    kval : Float
        Range of k-minus rates (per second) to use in simulation.
    time : Int
        Number of time steps to run the simulation.
    f_len : Int
        Filament length (subunits) to use in simulation.
    num_f : Int
        Formin number to use in simulation.
    k_plus : Float
        Rate of filament assembly (subunits/sec).
    TotalCables : Int
        Range of time to run the simulation.

    Returns
    -------
    taper_accum : DataFrame
        Returns an dataframe containing the maximum cable length at each time
        step and the parameters used during that simulation.

    '''
    #iniitalize an empty dataframe to store the results of the simulation
    taper_accum = pd.DataFrame()
    #iterate over the number of independent simulations to run
    for trajec in range(TotalCables):
        #compute the probability of filament removal for each
        #parameter set
        prob_var =  1 - np.exp(-1*(kval/(k_plus/f_len)))
        #initialize an empty dataframe        
        data = pd.DataFrame()        
        #Run the simluation and store the results.
        Lfila = FilaCompEvol(prob_var, Time, num_f)
        #Convert the arrays in the dataframe into columns where their position
        #is denoted 
        Lfila = pd.concat([Lfila, Lfila.pop("length").apply(pd.Series).add_prefix("")], axis=1)
        #replace the NaNs with zeros
        Lfila = Lfila.replace(np.nan, 0)
        #sum up the columns to get frequency of occupancy
        Lfila.loc['total'] = Lfila.sum()
        #transpose the dataframe and reset the index; need to remove the cell count #
        #as well
        data = Lfila.loc['total'][1:].T.reset_index()
        #add a column to record the length/position occupied
        data['dist'] = np.arange(0, data['total'].size)*(f_len*2.7/1000)
        data = data[['total', 'dist']]
        # dist_real = np.arange(0, 7.50595, 0.03869) #experimental distance interval
        # dist_real = pd.Series(dist_real, index=dist_real) #convert to a series for interpolation
        # data = data.values #convert to an array
        # data = pd.Series(data[:,0], index=data[:,1]) #convert to a series for interpolation
        # data = pd.DataFrame({'total':data, 'dist':dist_real}) #convert to DF
        # data = data.interpolate('index') #interpolate
        #record the parameters used in each simulation run        
        data['k_minus'] = kval * f_len
        data['iteration'] = trajec+1
        data['formin_num'] = num_f

        #data.drop_duplicates('dist', inplace=True)
        data.reset_index(inplace=True)
        taper_accum = taper_accum.append(data) #compile all of the data for each simulation
        #taper_accum = taper_accum.loc[taper_accum['dist'].isin(dist_real)]
    return taper_accum

#=============================================================================
#These functions are the analytic solutions to the theory for each cable
#behavior

#Cable trajectory=============================================================
import scipy.special as sc #import for function below
#below is a function to fit the expression for the cable trajectory
def analytic_traj(time, k_plus, k_minus, len_f, Nf):
    '''
    Parameters
    ----------
    time : Array, float
        Range of times to plot over.
    k_plus : Float
        Rate of filament assembly (per second).        
    k_minus : Float
        Rate of filament removal (per second).
    len_f : Float
        Length of filament (microns).
    Nf : Int
        Number of formins.

    Returns
    -------
    L_t : Array, float
        Cable length (microns) as a function of time (seconds).

    '''
    #compute the decay length, lambda
    lamb = (k_plus * len_f) / k_minus
    #compute the predicted cable length as a function of time
    L_t = lamb * (k_minus*time - 
                  np.exp((lamb/len_f)*Nf*np.exp(-1*k_minus*time))
                  * (sc.expi(-1*(lamb/len_f)*Nf)
                     - sc.expi(-1*(lamb/len_f)*Nf*np.exp(-1*k_minus*time))))
    return L_t

#Cable length distribution====================================================

#below is a function to fit the expression for the distribution for 
#cable lengths
def analytic_dist_cdf(L, Nf, lamb, len_f):
    '''
    Parameters
    ----------
    L : Array, float
        Range of lengths to compute the CDF/PDF.
    Nf : Int
        Number of formins.
    lamb : Float
        Decay length computed from k_plus, k_minus, and len_f.
    len_f : Float
        Length of filaments (microns).

    Returns
    -------
    prob_cdf : Array, float
        Array containing the CDF of a given length given the
        specified parameters.

    '''
    prob_L = pd.DataFrame() #initialize an empty dataframe
    prob_L['L'] = L #record the lengths used
    f = -1*(L/lamb) #compute outside of the function for ease of writing
    #Calculate the probability of each length
    prob_L['prob'] = Nf * np.exp(f - (lamb/len_f)*Nf*np.exp(f))
    #Compute the cumulative distributions for the analytic solution
    prob_L['cdf'] = np.cumsum(prob_L['prob'])
    prob_L['cdf_norm'] = prob_L['cdf'] / prob_L['cdf'].max()
    prob_cdf = prob_L['cdf_norm']

    return prob_cdf

def analytic_dist_pdf(L, Nf, lamb, len_f):
    '''

    Parameters
    ----------
    L : Array, float
        Range of lengths to compute the CDF/PDF.
    Nf : Int
        Number of formins.
    lamb : Float
        Decay length computed from k_plus, k_minus, and len_f.
    len_f : Float
        Length of filaments (microns).

    Returns
    -------
    prob_pdf : Array, float
        Array containing the PDF of a given length given the
        specified parameters.

    '''
    prob_L = pd.DataFrame() #initialize an empty dataframe
    prob_L['L'] = L #record the lengths used
    f = -1*(L/lamb) #compute outside of the function for ease of writing
    #Calculate the probability of each length
    prob_L['prob'] = Nf * np.exp(f - (lamb/len_f)*Nf*np.exp(f)) 
    #Compute the PDF for the analytic solution
    prob_pdf = prob_L['prob']

    return prob_pdf

#Cable density profile ======================================================= 

#below is a function to fit the expression for the cable tapering
def analytic_taper(L, F_0, k_minus_micron, k_plus_micron, len_f):
    '''
    Parameters
    ----------
    L : Array, float
        Range of lengths to compute the CDF/PDF.
    F_0 : Int
        Relative initial cable thickness (should be 1).
    k_minus_micron : Float
        Rate of filament removal (microns/sec).
    k_plus_micron : Float
        Rate of filament addition (microns/sec).        
    len_f : Float
        Length of filaments (microns).

    Returns
    -------
    F_l : Array, float
        Array containing the relative actin cable thickness as a function
        of length.
    '''
    
    F_l = F_0 * np.exp(-1 * (k_minus_micron*L)/(k_plus_micron*len_f))

    return F_l

#=============================================================================
#Define the parameters to be used during the simulations

#formin number to use in the simulations
Nf = 4 #number of formins
#set the range of filament lengths to simulate
f_len = 185 #subunits
#set the rate of assembly  
k_plus = 93 #subunits/sec
#set the k-minus rate
k_minus =  0.16 #per second
#set the total time to conduct the simulations
Time = 45 #steps
#set the number of trajectories to run for each simulation
TotalTraj = 500

#=============================================================================
#Indicate the parameters required for plotting the analytic solutions
time = np.linspace(0, 200, 201) #time, sec
L = np.linspace(0, 20, 100) #cable lengths, microns
k_plus_micron = 0.25 #um/sec
k_minus_sub = k_minus*f_len #subunits/sec
k_minus_micron = (k_minus_sub*2.7)/1e3 #um/sec
len_f = (f_len * (2.7/1e3)) #length of filaments in cable, microns
lamb = k_plus_micron / k_minus #expression for lambda, microns
#=============================================================================
#run the simulations for the specified parameters below
sim_data = pd.DataFrame(cable_sim(k_minus, Time, f_len, Nf, k_plus, TotalTraj))

taper_sim_data = pd.DataFrame(TaperSim(k_minus, Time, f_len, Nf, k_plus, TotalTraj))
#=============================================================================
#polish sim data for plotting

#find the final lengths from each trajectory
total_length = sim_data.groupby(['iteration','k_minus']).last().reset_index()
#compute the relative density of simulated cables
taper_sim_data['rel_total'] = taper_sim_data['total']/taper_sim_data['formin_num']
taper_sim_data['rel_norm'] = taper_sim_data['rel_total'] / taper_sim_data['rel_total'].max()
taper_sim_data = taper_sim_data.reset_index()

#=============================================================================
#Initiate plotting parameters
ft = 26
st = 'whitegrid'
cmap = 'viridis'
expt_color = '#C200FB'
sim_color = 'k'
theory_color = 'r'
#=============================================================================
#import the trajectory data to plot with simulation and theory results

#indicate the directory and file containing the data
data_traj = 'C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_trajectory_data\\'      
    
df_traj = pd.read_csv(data_traj + '200901_all_yBG12-Abp140Envy_trajectories.csv')

#group data by experimental replicate and time to calculate mean values
df_expt_mean = df_traj.groupby(['lifetime', 'n'],\
                                        sort=False).mean().reset_index()

#make an array of values to use in plotting the analytic result
traj_est_pars = [k_plus/f_len, k_minus, len_f, Nf]    

#=============================================================================
    
#plot the change in cable length over time for each replicate
with sns.axes_style(st):
    plt.figure(figsize=(6,6))
    #plot the simulation results            
    ax = sns.lineplot(x=sim_data['time'],
                      y=sim_data['length_rescale'],  
                      color= sim_color, lw=4, label='Simulation') 
    #plot the analytic result
    ax = sns.lineplot(x = time, y = analytic_traj(time, *traj_est_pars),
                      color=theory_color, lw=3, linestyle='--', label='Theory')     
    #plot the data
    ax = sns.lineplot(x=df_expt_mean['lifetime'],
                      y=df_expt_mean['neck_dist'],
                        color=expt_color, lw=4, label='Experiment')    
      
      
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10)) 
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    ax.tick_params(axis='y', colors='black',labelsize=ft)
    ax.tick_params(axis='x', labelsize=ft)
    plt.ylim([0, 10])
    plt.xlim([0, 45])        
    plt.legend(prop={'size': 20}) 
    plt.tight_layout()
    # plt.savefig('230706_sim_theory_traj.svg') 

#=============================================================================    
#import distribution data and add to the plots above
#import files to analyze
dist_data = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_length_variance_data\\"

#initalize data frame to append all data 
df_dist = pd.DataFrame()
#import data to dataframe
df_dist = pd.read_csv(dist_data + 
                  '210105_cdc28-13ts_t-8_t1_yBG12_yBG9_all_cable_analysis.csv')

df_dist = df_dist.loc[df_dist.strain=='yBG12', 'L'].values          #haploids


#make an array of values to use in plotting the analytic results
dist_pars = [Nf, lamb, len_f]   
#=============================================================================

#plot the cumulative distribution of cable length for each replicate

#combine the sim and expt data for plotting - this is a bit simpler
sim_dist = pd.DataFrame()
sim_dist['L'] =  total_length['length_rescale']
sim_dist['type'] = 'Sim'
expt_dist = pd.DataFrame()
expt_dist['L'] = df_dist
expt_dist['type'] = 'Expt'

frames = [expt_dist, sim_dist]

all_dist = pd.concat(frames)

cmap_small = [expt_color, sim_color]
with sns.axes_style(st):
    plt.figure(figsize=(6, 6))
    #Plot the simulation results   
    ax = sns.histplot(data = all_dist, x ='L', stat='probability', 
                      kde=True, hue='type',
                      fill=False,
                      palette=cmap_small, lw=0,
                      label='Simulation')
    

    ax.xaxis.set_major_locator(ticker.MultipleLocator(5)) 
    plt.xlabel(u'Cable length (${\mu}m$)',
                fontsize=ft)
    plt.ylabel('Probability', fontsize=ft) 
    plt.xlim([0,15])
    # plt.ylim([0,0.15])    
    ax.tick_params(axis='y', colors='black',labelsize=ft)
    ax.tick_params(axis='x', labelsize=ft)
    plt.legend([],[], frameon=False) 
    plt.tight_layout()
    # plt.savefig('231005_sim_expt_dist_wt_zoom.svg')
        
#=============================================================================
#Define the directory containing the anaylzed and compiled tapering data.
data_taper = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_intensity_profile_data\\final_data\\"
#Initialize an empty dataframe to write the data to.
df = pd.DataFrame()
#Read in the data to the empty dataframe.
df = pd.read_csv(data_taper + '230421_all_cable_traces_rescaled_trim.csv')
#Determine the length of each cable analyzed (i.e. the maximum length measured)
df['cable_length'] = df.groupby('cell')['Distance_(microns)'].transform('max')
#initialize dataframe to write cell-specific data to
df_cell = pd.DataFrame()
df_cell = df.groupby(['cell']).max().reset_index()
#Calculate the mean for each experimental replicate in each bin.
df_all_mean = df.groupby(['bin', 'n'], sort=True).mean().reset_index()
#make an array containing the variable for the analytic solution
taper_est_pars = [1, k_minus_micron, k_plus_micron, len_f]

#==============================================================================
#Plot the tapering data, analytic solution, and simulation results    

with sns.axes_style(st):
    plt.figure(figsize=(6, 6))
    #Plot the simulation results
    ax = sns.lineplot(x=taper_sim_data['dist'], y=taper_sim_data['rel_norm'],  
                      color=sim_color, lw=4, label='Simulation')

    #Plot the analytic result
    sns.lineplot(x = L, y = analytic_taper(L, *taper_est_pars), color=theory_color, \
                  lw=3, linestyle='--', label='Theory')
        
    #Plot the experimental data              
    ax = sns.lineplot(x= df['Distance_(microns)'].round(1),
                      y= df['int_norm_rescale'],lw=4, color=expt_color,
                      label='Experiment')

    ax.set_xlabel(u'Cable length (${\mu}m$)',
                fontsize=ft)
    ax.set_ylabel('Relative Cable Mass', fontsize=ft, color='black')
    plt.xlim([0, 6])
    ax.set_ylim([0, 1.1])
    ax.tick_params(axis='y', colors='black',labelsize=ft)
    ax.tick_params(axis='x', labelsize=ft)
    plt.legend(prop={'size': 20}) 
    
    plt.tight_layout()
    # plt.savefig('230706_sim_theory_taper.svg') 
    
#=============================================================================
#use the function below to compute the mean cable length for different
#combinations of parameters

def analytic_dist_pdf(L, Nf, lamb, len_f):
    '''

    Parameters
    ----------
    L : Array, float
        Range of lengths to compute the CDF/PDF.
    Nf : Int
        Number of formins.
    lamb : Float
        Decay length computed from k_plus, k_minus, and len_f.
    len_f : Float
        Length of filaments (microns).

    Returns
    -------
    prob_pdf : Array, float
        Array containing the PDF of a given length given the
        specified parameters.

    '''
    prob_L = pd.DataFrame() #initialize an empty dataframe
    prob_L['L'] = L #record the lengths used
    f = -1*(L/lamb) #compute outside of the function for ease of writing
    #Calculate the probability of each length
    prob_L['prob'] = Nf * np.exp(f - (lamb/len_f)*Nf*np.exp(f)) 
    #Compute the PDF for the analytic solution
    prob_pdf = prob_L['prob']

    return prob_pdf

#=============================================================================
#Compute the mean length (max P(L)) for each Nf and flen
#=============================================================================
#Indicate the parameters required for plotting the analytic solutions
L = np.linspace(0, 100, 1000) #cable lengths, microns
lamb = k_plus_micron / 0.15 #expression for lambda, microns  
Nf_sub = np.arange(2, 21, 2)
f_len_sub = np.arange(0.25, 2.6, 0.25)


df_Nf_flen = pd.DataFrame()

for flen in f_len_sub:
    for Nf in Nf_sub:
        df = pd.DataFrame()
        df['prob'] = analytic_dist_pdf(L, Nf, lamb, flen)
        df['L'] = L
        df['flen'] = flen*1000
        df['Nf'] = Nf
        df_Nf_flen = df_Nf_flen.append(df.loc[df['prob'].idxmax()],
                                       ignore_index=True) 
        
data = df_Nf_flen.pivot_table(index='Nf', columns='flen',
                            values='L').sort_values(by=['Nf'],
                                                           ascending=False)
                                                        
#=============================================================================
#Initiate plotting parameters
ft = 24
ft2 = 18
st = 'ticks'
cmap = 'viridis'
#=============================================================================
from matplotlib import ticker

#plot the results
with sns.axes_style(st):
    plt.figure(figsize=(8, 8))
    # sns.set_palette(cmap)
    ax = sns.heatmap(data, cmap="PRGn",
                      cbar_kws={'label': u'Mean cable length (${\mu}m$)'},
                      center=4.3,
                      # vmin = 0,
                      # vmax = 8.6,
                      annot=True, fmt=".1f",
                      yticklabels=1, linewidths=1, linecolor='k')    
            
    ax.legend([],[], frameon=False)
    cbar = ax.collections[0].colorbar
    cbar.ax.yaxis.label.set_size(ft2)
    cbar.ax.tick_params(labelsize=ft2)
    
    plt.ylabel(u'Formin number', fontsize=ft)
    # plt.ylabel(u'Filament length (nm)', fontsize=ft)
    
    plt.xlabel(u'Filament length (nm)', fontsize=ft)
    # plt.xlabel(u'Formin number', fontsize=ft)    
    plt.xticks(rotation=45, fontsize=ft)
    plt.yticks(rotation=0, fontsize=ft)
    plt.tight_layout(h_pad=1, w_pad=1)
    # plt.title('Tapering/Trajectory BIC', fontsize=ft)
    # plt.savefig('231009_theory_Nf_Lf_ratio.svg')   
    
    