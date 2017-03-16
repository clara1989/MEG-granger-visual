# -*- coding: utf-8 -*-
"""
Created on Mon Jan 09 09:25:11 2017

@author: seymourr
"""

###############################################################################
#
# V1 --> V4 Granger Causality
#
###############################################################################

cd C:\Users\seymourr\Google Drive\RS PhD Documents\Results\Granger

import csv
import numpy as np
from matplotlib import pyplot as plt

ff_ASD_file = 'D:/ASD_Data/Group/GC/ASD_V1_V4_ff.csv';
fb_ASD_file = 'D:/ASD_Data/Group/GC/ASD_V1_V4_fb.csv';

## Extract the important information
def extract_granger(filename):
    temp_dict = []        
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            temp_dict.append(row)
            
    temp_dict = np.array(temp_dict).astype(np.float)
    return temp_dict   
  
ff_ASD = extract_granger(ff_ASD_file) 
fb_ASD = extract_granger(fb_ASD_file)  

# Calculate Means
ff_ASD_mean = ff_ASD.mean(axis = 0)
fb_ASD_mean = fb_ASD.mean(axis = 0)

## Calculate CIs
import scipy as sp
import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return h, m-h, m+h
    
CI_fb_ASD = mean_confidence_interval(fb_ASD,confidence=0.95)
CI_ff_ASD = mean_confidence_interval(ff_ASD,confidence=0.95)


## Plot
x = range(1,141,1)

import seaborn as sns
sns.set(style="white")


plt.figure()
plt.plot(x, fb_ASD_mean, 'k-',color='#1B2ACC')
plt.plot(x, ff_ASD_mean, 'k-',color='#D13D1D')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.fill_between(x,fb_ASD_mean+CI_fb_ASD[0],fb_ASD_mean-CI_fb_ASD[0],color='#b9cfe7', edgecolor='')
plt.fill_between(x,ff_ASD_mean+CI_ff_ASD[0],ff_ASD_mean-CI_ff_ASD[0],color='#ffa5a5', edgecolor='')
plt.xlabel('Frequency (Hz)',fontsize=20)
plt.ylabel('Granger Causality Value',fontsize=20)
plt.legend(['Feedback','Feedforward'],fontsize=20)
plt.savefig("ASD_GC.png", dpi=600)
plt.show()

###############################################################################
#
# DAI
#
###############################################################################

DAI_file_pilot = 'D:/ASD_Data/Group/GC/DAI_pilot.csv';
DAI_file_ASD = 'D:/ASD_Data/Group/GC/DAI_ASD.csv';
  
DAI_pilot = extract_granger(DAI_file_pilot) 
DAI_ASD = extract_granger(DAI_file_ASD) 

DAI_pilot_mean = DAI_pilot.mean(axis = 0)
DAI_ASD_mean = DAI_ASD.mean(axis = 0) 

CI_DAI_pilot = mean_confidence_interval(DAI_pilot,confidence=0.95)
CI_DAI_ASD = mean_confidence_interval(DAI_ASD,confidence=0.95)

x = range(1,141,1)
zero_line = [0] * 140 

import seaborn as sns
sns.set(style="white")

plt.figure()
plt.plot(x, DAI_pilot_mean, 'k-',color='#009e02')
plt.plot(x, DAI_ASD_mean, 'k-',color='#f95800')
plt.plot(x, zero_line, 'k-',color='#000000')
plt.fill_between(x,DAI_pilot_mean+CI_DAI_pilot[0],DAI_pilot_mean-CI_DAI_pilot[0],color='#A5DB94', edgecolor='')
plt.fill_between(x,DAI_ASD_mean+CI_DAI_ASD[0],DAI_ASD_mean-CI_DAI_ASD[0],color='#ffdecc', edgecolor='')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Granger Causality Value')
plt.legend(['Pilot DAI','Control DAI'])
plt.show()

