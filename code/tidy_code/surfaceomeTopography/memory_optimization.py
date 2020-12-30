#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:08:56 2019

@author: southk

memory_optimization.py

functions for optimizing memory usage when reading a large file into pandas

adapted from https://www.dataquest.io/blog/pandas-big-data/

"""
import pandas as pd

def mem_usage(pandas_obj):
    
    if isinstance(pandas_obj,pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
        
    else: # we assume if not a df it's a series
        usage_b = pandas_obj.memory_usage(deep=True)
        
    usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
    return "{:03.2f} MB".format(usage_mb) 


def optimize_int(full):
    
    full_int = full.select_dtypes(include=['int'])
    converted_int = full_int.apply(pd.to_numeric,downcast='unsigned')
    return converted_int

    
def optimize_obj(full):
    
    full_obj = full.select_dtypes(include=['object']).copy()
    
    converted_obj = pd.DataFrame()

    for col in full_obj.columns:
        num_unique_values = len(full_obj[col].unique())
        num_total_values = len(full_obj[col])
        if num_unique_values / num_total_values < 0.5:
            converted_obj.loc[:,col] = full_obj[col].astype('category')
        else:
            converted_obj.loc[:,col] = full_obj[col]
    
    return converted_obj
            
def optimize_all(full):
    optimized_full = full.copy()
    
    converted_int = optimize_int(full)
    optimized_full[converted_int.columns] = converted_int
    
    converted_obj = optimize_obj(full)
    optimized_full[converted_obj.columns] = converted_obj
    
    return optimized_full
    
    
    
    
