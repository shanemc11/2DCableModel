# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 12:45:44 2023

@author: shane
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import random as rand
from statsmodels.distributions.empirical_distribution import ECDF
import seaborn as sns 
sns.set()
#import curve_fit from scipy
from scipy.optimize import curve_fit

#=============================================================================
#Goal: This script contains the functions required to simulate cable
#assembly using the cable treadmilling model. It will iterate through the 
#specified ranges of parameter values and output these results as csv files
#to the specified directory. Due to the large range of parameters, this script
#requires a few days to run.
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
        all_filaments = pd.concat([all_filaments, Lfilaments], axis=0,
                                  ignore_index=True)
        
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
        prob_var =   1 - np.exp(-1*(kval/k_plus))
        #initialize an empty dataframe
        data = pd.DataFrame()
        #run simulation        
        Lcable = CompCableEvol(prob_var, time, num_f, f_len)
        #convert array values to lengths (microns)        
        data['length'] = Lcable * (f_len*2.7/1000)
        #rescale lengths to start growth at zero length
        data['length_rescale'] = (data['length']-data['length'][0])
        #scale the time into seconds
        data['time'] = np.arange(0, time, 1)*(f_len)/(k_plus*f_len)
        #record the parameters used in each simulation run        
        data['k_minus'] = kval * f_len
        data['iteration'] = i+1
        #append the data for each iteration
        data_accum = pd.concat([data_accum, data], axis=0,
                               ignore_index=True)
        
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
        prob_var =  1 - np.exp(-1*(kval/k_plus))
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
        # dist_real = np.arange(0, 10.75, 0.02136) #experimental distance interval
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
        taper_accum = pd.concat([taper_accum, data], axis=0,
                               ignore_index=True)
        #taper_accum = taper_accum.loc[taper_accum['dist'].isin(dist_real)]
    return taper_accum
#=============================================================================
#define the parameters to be used during the simulations

#range of formin number to use in the simulations
NoCables = 4 #number of formins
#set the range of filament lengths to simulate
f_len = 185 #subunits
#set the rate of assembly  
k_plus = 0.50 #subunits/sec
# indicate the range of k-minus values to expolore 
k_minus =  0.16 #per second
#set the total time to conduct the simulations
Time = 120 #time steps
#set the number of trajectories to run for each simulation
TotalTraj = 250
#indicate the ratio of cell lengths
ratio = 1.9
#=============================================================================
#Indicate the parameters required for plotting the analytic solutions
# time = np.linspace(0, 200, 201) #time, sec
# L = np.linspace(0, 20, 100) #cable lengths, microns
# k_plus_micron = 0.25 #um/sec
# k_minus_sub = k_minus*f_len #subunits/sec
# k_minus_micron = (k_minus_sub*2.7)/1e3 #um/sec
# len_f = (f_len * (2.7/1e3)) #length of filaments in cable, microns
# lamb = k_plus_micron / k_minus #expression for lambda, microns
#=============================================================================

#run the simulations given the parameters specified above for small cells
sim_small = pd.DataFrame(cable_sim(k_minus, Time, f_len, NoCables, k_plus, TotalTraj))
#round the k_minus values
sim_small['k_minus'] = sim_small['k_minus'].round(0)

taper_sim_small = pd.DataFrame(TaperSim(k_minus, Time, f_len, NoCables, k_plus, TotalTraj))
#round the k_minus values
taper_sim_small['k_minus'] = taper_sim_small['k_minus'].round(0)

#=============================================================================
#polish sim data for plotting

#find the final lengths from each trajectory
total_length_small = sim_small.groupby(['iteration','k_minus']).last().reset_index()
#compute the relative density of simulated cables
taper_sim_small['rel_total'] = taper_sim_small['total']/taper_sim_small['formin_num']
taper_sim_small['rel_norm'] = taper_sim_small['rel_total'] / taper_sim_small['rel_total'].max()
taper_sim_small = taper_sim_small.reset_index()
#drop lengths that are larger than the longest measured cable
taper_sim_small = taper_sim_small[taper_sim_small['dist'] < 5.75] 

#=============================================================================
#This section test the dependence of f_len on cell length

#run the simulations given the parameters specified above for large cells
sim_large = pd.DataFrame(cable_sim(k_minus/ratio, Time, f_len*ratio,
                                    NoCables, k_plus/ratio, TotalTraj))
#round the k_minus values
sim_large['k_minus'] = sim_large['k_minus'].round(0)

taper_sim_large = pd.DataFrame(TaperSim(k_minus/ratio, Time, f_len*ratio,
                                        NoCables, k_plus/ratio, TotalTraj))
#round the k_minus values
taper_sim_large['k_minus'] = taper_sim_large['k_minus'].round(0)

#=============================================================================
#This section test the dependence of k-minus on cell length
#run the simulations given the parameters specified above for large cells
# sim_large = pd.DataFrame(cable_sim(k_minus/ratio, Time, f_len,
#                                     NoCables, k_plus, TotalTraj))
# #round the k_minus values
# sim_large['k_minus'] = sim_large['k_minus'].round(0)

# taper_sim_large = pd.DataFrame(TaperSim(k_minus/ratio, Time, f_len,
#                                         NoCables, k_plus, TotalTraj))
# #round the k_minus values
# taper_sim_large['k_minus'] = taper_sim_large['k_minus'].round(0)

#=============================================================================
#polish sim data for plotting

#find the final lengths from each trajectory
total_length_large = sim_large.groupby(['iteration','k_minus']).last().reset_index()
#compute the relative density of simulated cables
taper_sim_large['rel_total'] = taper_sim_large['total']/taper_sim_large['formin_num']
taper_sim_large['rel_norm'] = taper_sim_large['rel_total'] / taper_sim_large['rel_total'].max()
taper_sim_large = taper_sim_large.reset_index()

#=============================================================================
#Initiate plotting parameters
ft = 26
st = 'whitegrid'
cmap = ['#6e7f80', '#000000']
expt_color = '#C200FB'
sim_color = 'k'
theory_color = 'r'
#=============================================================================
#import the trajectory data to plot with simulation and theory results

#setup the data directory
datadir = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_trajectory_data\\"

#initalize data frame to append all data 
df_t0 = pd.DataFrame()

df_t8 = pd.DataFrame()

all_df = pd.DataFrame()


#read in the summary data files to compare t0 to t8 cells
df_t0 = pd.read_csv(datadir + \
                    "200826_t0_all_cable_extension_rate_analysis.csv")
df_t0 = df_t0.drop(columns=['cell', 'time'])

df_t8 = pd.read_csv(datadir + \
                    "200826_t8_all_cable_extension_rate_analysis.csv")
df_t8 = df_t8.drop(columns=['cell', 'time'])

    
#calculate means for each timepoint for each replicate
df_t0_expt_mean = df_t0.groupby(['lifetime', 'n'], \
                                sort=False).mean().reset_index()
    
df_t8_expt_mean = df_t8.groupby(['lifetime', 'n'],\
                                       sort=False).mean().reset_index()    

#make an array of values to use in plotting the analytic result
# traj_est_pars = [k_plus/f_len, k_minus, len_f, Nf]    

#=============================================================================
    
#plot the change in cable length over time for each replicate
with sns.axes_style(st):
    plt.figure(figsize=(6,6))
    #plot the simulation results            
    ax = sns.lineplot(x=sim_small['time'],
                      y=sim_small['length_rescale'],  
                      color='#000000', lw=4, label='Small')
    #plot the simulation results            
    ax = sns.lineplot(x=sim_large['time'],
                      y=sim_large['length_rescale'],  
                      color='#6e7f80', lw=4, label='Large')    
    #plot the analytic result
    # ax = sns.lineplot(time, analytic_traj(time, *traj_est_pars),
    #                   color=theory_color, lw=3, linestyle='--', label='Theory')     
    #plot the data
    #plot the mean and 95%CI for the replicates   
    ax = sns.lineplot(x=df_t0_expt_mean['lifetime'], \
                      y=df_t0_expt_mean['neck_dist'], \
                          errorbar=('ci', 95), label='cdc28-13, uninduced',\
                              color='#7fb800', lw=3)
    
    ax = sns.lineplot(x=df_t8_expt_mean['lifetime'], \
                      y=df_t8_expt_mean['neck_dist'],\
                      color='#ffb400', errorbar=('ci', 95),\
                          label = 'cdc28-13, induced',\
                          lw=3)    
      
      
    plt.xlabel('Extension time (sec)', fontsize=ft)    
    plt.ylabel(u'Cable length (${\mu}m$)', fontsize=ft)
    plt.rc('xtick', labelsize=ft)
    plt.rc('ytick', labelsize=ft)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(15))     
    ax.tick_params(axis='y', colors='black',labelsize=ft)
    ax.tick_params(axis='x', labelsize=ft)
    plt.ylim([0, 10])
    plt.xlim([0, 60])        
    # plt.legend(prop={'size': 10})
    plt.legend([],[], frameon=False) 
    plt.tight_layout()
    # plt.savefig('230706_sim_expt_traj_cdc28_predictions_kminus_scaling.svg')

#=============================================================================    
#import distribution data and add to the plots above
data_dist = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_length_variance_data\\"

#initalize data frame to append all data 
df_dist = pd.DataFrame()
#import data to dataframe
df_dist = pd.read_csv(data_dist + 
                  '210105_cdc28-13ts_t-8_t1_yBG12_yBG9_all_cable_analysis.csv')   

t0_dist = df_dist[(df_dist['strain'] == 'cdc28-13ts, t0')].reset_index() 
t8_dist = df_dist[(df_dist['strain'] == 'cdc28-13ts, t8')].reset_index()  
 
frames = [t0_dist, t8_dist]

dist_all = pd.concat(frames).reset_index()

#make an array of values to use in plotting the analytic results
# dist_pars = [Nf, lamb, len_f]   
#=============================================================================

#combine the sim and expt data for plotting - this is a bit simpler
sim_dist_t0 = pd.DataFrame()
sim_dist_t0['L'] =  total_length_small['length_rescale']
sim_dist_t0['type'] = 'Sim'
expt_dist_t0 = pd.DataFrame()
expt_dist_t0['L'] = t0_dist['L']
expt_dist_t0['type'] = 'Expt'

frames = [expt_dist_t0, sim_dist_t0]

t0_all_dist = pd.concat(frames)
cmap_small = ['#7fb800', 'black']
with sns.axes_style(st):
    plt.figure(figsize=(6, 6))
    #Plot the simulation results   
    ax = sns.histplot(data = t0_all_dist, x ='L', stat='probability', 
                      kde=True, hue='type',
                      fill=False,
                      palette=cmap_small, lw=0,
                      label='Simulation')
    

    ax.xaxis.set_major_locator(ticker.MultipleLocator(5)) 
    plt.xlabel(u'Cable length (${\mu}m$)',
                fontsize=ft)
    plt.ylabel('Probability', fontsize=ft) 
    plt.xlim([0,20])
    # plt.ylim([0,0.15])    
    ax.tick_params(axis='y', colors='black',labelsize=ft)
    ax.tick_params(axis='x', labelsize=ft)
    plt.legend([],[], frameon=False) 
    plt.tight_layout()
    # plt.savefig('231004_sim_expt_dist_cdc28t0_predictions_flen.svg')
    
#==============================================================================    
sim_dist_t8 = pd.DataFrame()
sim_dist_t8['L'] =  total_length_large['length_rescale']
sim_dist_t8['type'] = 'Sim'
expt_dist_t8 = pd.DataFrame()
expt_dist_t8['L'] = t8_dist['L']
expt_dist_t8['type'] = 'Expt'

frames = [expt_dist_t8, sim_dist_t8]

t8_all_dist = pd.concat(frames)
cmap_large = ['#ffb400', 'black']
with sns.axes_style(st):
    plt.figure(figsize=(6, 6))
    #Plot the simulation results
    ax = sns.histplot(data = t8_all_dist, x ='L', stat='probability', 
                      kde=True, hue='type',
                      fill=False,
                      palette=cmap_large, lw=0)
    

    ax.xaxis.set_major_locator(ticker.MultipleLocator(5)) 
    plt.xlabel(u'Cable length (${\mu}m$)',
                fontsize=ft)
    plt.ylabel('Probability', fontsize=ft) 
    plt.xlim([0,30])
    # plt.ylim([0,0.15])    
    ax.tick_params(axis='y', colors='black',labelsize=ft)
    ax.tick_params(axis='x', labelsize=ft)
    plt.legend([],[], frameon=False) 
    plt.tight_layout()
    # plt.savefig('231004_sim_expt_dist_cdc28t8_predictions_flen.svg')
            
#=============================================================================
#Define the directory containing the anaylzed and compiled tapering data.
taper_dir = "C:\\Users\\shane\\OneDrive\\Desktop\\temp_src\\cable_geometry_paper_data\\cable_intensity_profile_data\\final_data\\"
#Initialize an empty dataframe to write the data to.
df = pd.DataFrame()
#Read in the data to the empty dataframe.
df = pd.read_csv(taper_dir + '230510_cdc28_all_cable_traces_rescaled_clean.csv')
# df = df.drop(columns=['cell'])

#Determine the length of each cable analyzed (i.e. the maximum length measured)
df['cable_length'] = df.groupby('cell')['Distance_(microns)'].transform('max')
#initialize dataframe to write cell-specific data to
df_cell = pd.DataFrame()
df_cell = df.groupby(['cell']).max().reset_index()
#Calculate the mean for each experimental replicate in each bin.
df = df.drop(columns=['cell'])

df_all_mean = df.groupby(['bin', 'n'], sort=True).mean().reset_index()
#make an array containing the variable for the analytic solution
# taper_est_pars = [1, k_minus_micron, k_plus_micron, len_f]

#group all of the data from bin 1 into a single dataframe
bin1 = df[(df['bin'] == 1)].reset_index()  
bin1 = bin1.groupby(['n', 'Distance_(microns)']).mean().reset_index()
#group all of the data from bin 3 into a single dataframe
bin3 = df[(df['bin'] == 3)].reset_index()  
bin3 = bin3.groupby(['n', 'Distance_(microns)']).mean().reset_index()


#==============================================================================
#Plot the tapering data, analytic solution, and simulation results    

# #parse the simulations based on k-minus rate
# taper_small = taper_sim_data[(taper_sim_data['k_minus'] == 
#                               (f_len*k_minus_var[1]).round(0))]
# #drop lengths that are larger than the longest measured cable
# taper_small = taper_small[taper_small['dist'] < 5.75] 

# taper_long = taper_sim_data[(taper_sim_data['k_minus'] == 
#                              (f_len*k_minus_var[0]).round(0))]


with sns.axes_style(st):
    plt.figure(figsize=(6, 6))
    #Plot the simulation results
    ax = sns.lineplot(x=taper_sim_small['dist'], y=taper_sim_small['rel_norm'],
                      color=cmap[1], lw=4, label='Simulation')
    
    ax = sns.lineplot(x=taper_sim_large['dist'], y=taper_sim_large['rel_norm'],
                      color=cmap[0], lw=4, label='Simulation')
    

    #Plot the analytic result
    # sns.lineplot(L, analytic_taper(L, *taper_est_pars), color=theory_color, \
    #              lw=3, linestyle='--', label='Theory')
        
    #Plot the experimental data              
    ax = sns.lineplot(x= bin1['Distance_(microns)'].round(1),
                      y= bin1['int_norm_rescale'], lw=4, color='#7fb800',
                      label='Cdc28-13ts, t=0')
    
    ax = sns.lineplot(x= bin3['Distance_(microns)'].round(1), 
                      y= bin3['int_norm_rescale'], lw=4, color='#ffb400',
                      label='Cdc28-13ts, t=8')        


    ax.set_xlabel(u'Cable length (${\mu}m$)',
                fontsize=ft)
    ax.set_ylabel('Relative Cable Mass', fontsize=ft, color='black')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2)) 
    plt.xlim([0, 10])
    ax.set_ylim([0, 1.1])
    ax.tick_params(axis='y', colors='black',labelsize=ft)
    ax.tick_params(axis='x', labelsize=ft)
    # plt.legend(prop={'size': 10}) 
    plt.legend([],[], frameon=False)    
    plt.tight_layout()
    # plt.savefig('230706_sim_expt_taper_cdc28_predictions_kminus_scaling.svg') 
    
  