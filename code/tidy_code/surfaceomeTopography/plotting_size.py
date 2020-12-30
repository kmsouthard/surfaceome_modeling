#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:37:04 2019

size distribution ploting functions

@author: southk
"""

from pylab import *
from scipy.optimize import curve_fit



def gauss(x,mu,sigma,A):
    return A*exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def trimodal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)

def quadramodal(x,mu1,sigma1,A1,mu2,sigma2,A2,mu3,sigma3,A3, mu4,sigma4,A4):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)+gauss(x,mu3,sigma3,A3)+gauss(x,mu4,sigma4,A4)


def gauss_hist_fit(df, fit = quadramodal, nbins=44, max_height=220, peaks =[40, 70, 125, 200]):
    
    data=df['max_dim'].tolist()
    total = df['max_dim'].count()
    
    peak_hieghts = [total*0.125, total*0.050, total*0.025, total*0.0125]
    
    y,x,_=hist(data,nbins,alpha=.3,label='data', range = (0,max_height))
    
    x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
    
    expected=(peaks[0],peak_hieghts[0],20,peaks[1],peak_hieghts[1],20,
              peaks[2],peak_hieghts[2],20,peaks[3],peak_hieghts[3],10)
    params,cov=curve_fit(fit,x,y,expected)
    sigma=sqrt(diag(cov))
    plot(x,fit(x,*params),color='red',lw=3,label='gaussian')
    legend()
    print(params,'\n',sigma)
    
    
    
def multi_hist(dfs, labels, col = 'max_dim', max_height = 220, normalize = True):
    
    plt.figure();
    
    count = 0
    
    for label in labels:
        dfs[count].name = label
        count += 1
    
    for df in dfs:
        
        df[col].plot.hist(alpha=0.5,
                           bins=range(0, max_height, 5), density = normalize, label = df.name)
        
        legend()
        
def multi_hist_cumulative(dfs, labels, col = 'max_dim', max_height = 220):
    
    plt.figure();
    
    count = 0
    
    for label in labels:
        dfs[count].name = label
        count += 1
    
    
    for df in dfs:

        df[col].plot.hist(alpha=0.5,
                           bins=range(0, max_height, 1), cumulative = True, histtype = 'step', label = df.name)