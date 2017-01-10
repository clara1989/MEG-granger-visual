# -*- coding: utf-8 -*-
"""
Created on Mon Jan 09 09:25:11 2017

@author: seymourr
"""

import csv
import numpy as np
from matplotlib import pyplot as plt

ff_ASD_file = 'D:/ASD_Data/Group/GC/ASD_ff.csv';
fb_ASD_file = 'D:/ASD_Data/Group/GC/ASD_fb.csv';

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
sns.set(style="whitegrid")

plt.figure()
plt.plot(x, fb_ASD_mean, 'k-',color='#1B2ACC')
plt.plot(x, ff_ASD_mean, 'k-',color='#D13D1D')
plt.fill_between(x,fb_ASD_mean+CI_fb_ASD[0],fb_ASD_mean-CI_fb_ASD[0],color='#b9cfe7', edgecolor='')
plt.fill_between(x,ff_ASD_mean+CI_ff_ASD[0],ff_ASD_mean-CI_ff_ASD[0],color='#ffa5a5', edgecolor='')
plt.xlabel('Frequency(Hz)')
plt.ylabel('GC Value')
plt.legend(['Feedback','Feedforward'])
plt.show()


