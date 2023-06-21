#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 15:32:14 2020

@author: benjamincolmey
"""
#%% import needed packages, 

# defining values, always keep at the top
import matplotlib.pyplot as plt
import numpy as np
import random 
import fast_histogram
from fast_histogram import histogram1d
from scipy.stats import poisson



Nsamp = int(1e6);# refers to the total number of pulses that will be produced. It must be an interger
PulsesPerSecond = int(76e6) # gives the number of pulses per second, or frequency, again as an integer
ElapsedTime = (Nsamp)/(PulsesPerSecond) # This is total time since laser is on in seconds, should be equal tp largest value in DectorTimes
Efficiency = 100 #is a global efficiency, which will reduce the number of photons to reach either detector by a percentage (must be written as an integer from 1-100)
Tau = 1 ## time for exponential to decay to 1/e times its initial value
BeamsplitterRatio = 0.5 #determines the given reflectance vs transmittance of our beamsplitter, written here as a decimal from 0-1
time_bw_pulses = 13.0 #(ns) is the time between successive excitation pulses in nanoseconds
peakCenter = 0
lamb = 1 #parameter which we will use as part of a Poissonian distribution model
sigma = 0.001 ## for Gaussian
mean = 100 ## for Gaussian

deadtime = 39 #ns
Proximity = 4.5  #(ns) how close to t= 0 is considered a second photon emission?
DarkCountRate = 200 # per second
DarkConst = DarkCountRate/PulsesPerSecond


Detect = True
Gauss = False
QuantumDot =  True
Dead = False
MultiEmission = False
Plot = True
AutoCorrelation = True
DarkCounts = False
TriggerPulseWidth = (1.5e-2) #150 fs
ProbOfSecondEmission  = 1 #on the order of 1 for approx 50/50 second emission once first emission falls before end of trigger pulse

#%% 
##This section creates a stream of photons meant to simulates our quantum dots, there are several parameters which can be turned
##on or off such as Guassian or Poissonian jitter, deadtime, efficiency, beamsplitter ratio

if QuantumDot:
    
    time = []
    TotalTime= []
    PhotonNumb = range(Nsamp)
    Detector1time = []
    Detector2time = []

     #creating empty arrays to be filled later
   
    #creating stream of photons
    for i in PhotonNumb:
        #creating the cumulative distribution function
        ro = random.random()
        r1 = random.random()
        time.append(-(1/Tau)*np.log(1 - ro))
        TotalTime.append((i+1)*13)
        
        if (DarkCounts == True):
            
            Dark = np.random.poisson(DarkConst,1)[0]  
           # d = 0
            if (Dark == 1):
                print('1 Dark Count')
                #d += 1
                shift = (random.random()-0.5)*13
                r6 = random.random()
            #print('Dark Counts =',  d)
                
                if (r6 < BeamsplitterRatio):
                #working in the assymmetric beamsplitter            
                    #detector 1 Arrivals corresponds to a value of 1
                    Detector1time.append(time[i] + TotalTime[i-1]+shift)
                    Detector2time.append(0)
                    
                    #keeping track of detetction time, which is total time + emission time for each photon 
                else:
                    #detector 1 Arrivals corresponds to a value of 1
                    Detector2time.append(time[i] + TotalTime[i-1]+shift)
                    Detector1time.append(0)
                    
        if Gauss:
                
            r2=random.random()
            r3=random.random()
            time.append(time[i]+ sigma*np.sqrt(2*-np.log(1-r2))*np.cos(2*np.pi*r3))
      
        if (r1 < BeamsplitterRatio):
            Detector1time.append(time[i] + TotalTime[i])
            Detector2time.append(0)
        else:
            Detector2time.append(time[i] + TotalTime[i])
            Detector1time.append(0)

        if MultiEmission:
            
            if (time[i] <= TriggerPulseWidth):
                r4 = random.random()
                r5 = random.random()
               
                if ((TriggerPulseWidth-time[i])*ProbOfSecondEmission*20000 > r4):
                    #print('Double detection event')
                    r = random.random()
                    time1 = (-(1/Tau)*np.log(1 - r)) #second photon lifetime 
                        
                    if (r5 < BeamsplitterRatio):
                        Detector1time.append(time[i] + time1 + TotalTime[i-1])
                        Detector2time.append(0)
                    else:                    
                        Detector2time.append(time[i] + time1 + TotalTime[i-1])
                        Detector1time.append(0)


#%% Processing dead times
if Dead:
    def detector_dead(times, deadtime):
        # if list is empty, doing nothing
        if len(times) == 0:
            return
        # otherwise adding first element with deadtime to get the time in which system will
        # be free to detect other values again
        current = deadtime + times[0]
        # looping through each of the remaining values
        for i in range(1, len(times)):
            # if data[i]<current, it means the deadtime for previous value is not completed
            # so we set element at index i to 0
            if times[i] < current:
                times[i] = 0
            # otherwise, it means previous deadtime time is completed, so we add deadtime
            # value to element at index i to get the next time the system will be free again.
            else:
                current = times[i] + deadtime
    detector_dead(Detector2time, deadtime) 
    detector_dead(Detector1time, deadtime)  
      
#%%Counting coincidence events

if Detect:
    #options:
    NewNSamp = len(Detector1time)
    numNearestEvents = 10; #number of nearest events to compare in each direction; 10 means 10 to left, 10 more to right.
    percent_updates = 5#percent; display progress update every X percent in console.
    update_every_X = int(NewNSamp*(percent_updates/100));#convert % to number of entries per X percent!
    mute_prints =  True;
    #storage of delays:
    Detector1delays = []; #for storing delays calculated using detector 1 first
    
    for i in range(NewNSamp):
        
        if (Detector1time[i] != 0): 
          
            #IF detector 1 Arrivals, look at <=21 nearest possible times for det2.
            
            #within 10 elements of the list START, fix left 'fencepost' at start.
            if i < numNearestEvents:
                Detector2time_slice = Detector2time[:i+numNearestEvents+1];
            #within 10 elements of the list END, fix left 'fencepost' at end of list.
            elif i > (NewNSamp-numNearestEvents):
                Detector2time_slice = Detector2time[i-numNearestEvents:];
            else: #comfortably away from the list edges
                Detector2time_slice = Detector2time[i-numNearestEvents:i+numNearestEvents+1];
            
            ### trimming & storage of delay times, common to all 3 above statements in if-elif-else
            if not(mute_prints): #print to check comparing to 2*N+1 closest possible values
                print("len det2timeslice pretrim: {0}".format(len(Detector2time_slice)))
            
            #trim zeros
            Detector2time_slice = [det2time for det2time in Detector2time_slice if det2time!=0];
            if not(mute_prints): #print to see pulled times
                print(i, np.round(Detector2time_slice,2))
            if (Detector1time[i] != 0):
                #calculate and store delay
                delay_times_to_store = np.array(Detector2time_slice) - Detector1time[i]
                if not(mute_prints): #print to see stored delays
                    print(i, delay_times_to_store)
                Detector1delays.append(delay_times_to_store)
                
        #display percent completed!
        if i%update_every_X == 0:
            percent_through = i/NewNSamp * 100;
            print("{0:3.0f}% complete processing delays for histogram.".format(percent_through))
    print("100% complete processing delays for histogram.".format(percent_through))
            
    
    #function to flatten a list of any kind of iterable -- https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    flatten = lambda l: [item for sublist in l for item in sublist]
    Detector1delays_flat = flatten(Detector1delays)
    print("\n{0} delays stored for histogram.".format(len(Detector1delays_flat)))
    
#%%Calculating g2

if AutoCorrelation:  
    
    ZeroCounts = 1
    #Creating an empty integer variable that will be added to
    for i in range(len(Detector1delays_flat)):
        if (Detector1delays_flat[i] < Proximity and Detector1delays_flat[i] > -Proximity):

        # looping through all delays that fall below or above the proximity time chosen above  
            ZeroCounts += 1
            #print('1 Zero Count Found')
    
    NonZeroCounts = len(Detector1delays_flat)-ZeroCounts 
    #total counts that dont fall at timedelay = 0
    
    PeakNumber = numNearestEvents * 2
    AveCountsPerPeak = NonZeroCounts/PeakNumber
    
    g_two =  ZeroCounts/AveCountsPerPeak
    g_two = round(g_two, 3)
    print('Zero Counts=', ZeroCounts)
    print('Non Zero Counts = ', NonZeroCounts)
    print('Average Counts Per Peak =', AveCountsPerPeak)
    print('g^2(0) = ', g_two)
     




#%% Plot histogram 

if Plot: 
    BinsPerInteger = 0.5
    Bins = np.arange(round(min(Detector1delays_flat)) - 0.5, round(max(Detector1delays_flat)) + 0.5, BinsPerInteger)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Generate the histogram and plot it
    n, bins, patches = ax.hist(Detector1delays_flat, bins=Bins, color="b", alpha=0.5)
    
    ax.text(0.08, 0.98, f"Total Detection events: {n.sum()}", horizontalalignment='left', verticalalignment='top', transform=ax.transAxes,fontsize=16)
    ax.text(0.08, 0.93, f"g^2(0) =  {g_two}", horizontalalignment='left', verticalalignment='top', transform=ax.transAxes,fontsize=16)
    ax.xaxis.set_major_locator(plt.MultipleLocator(int(time_bw_pulses)))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(Proximity))
    
    ax.set_xlabel('Time Delay (ns)', fontsize=15)
    ax.set_ylabel('Coincidence Events', fontsize=15)
    ax.set_xlim([5 * -time_bw_pulses, 5 * time_bw_pulses])
    ax.set_ylim(0, None)

    # Display the plot
    plt.tight_layout()
    plt.show()

