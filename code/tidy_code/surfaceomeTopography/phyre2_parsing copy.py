#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 14:45:17 2019

@author: southk

Processing Phyre2 output Data:

combining protein measurements and Phyre2 results annotation
"""

import pandas as pd
import numpy as np
import os
import re



def read_summaryinfo(path, pattern):
    #read in summaryinfo or summaryinfo_extended files from phyre2 results
    all_files = [(x+'/'+pattern) for x in os.listdir(path) if os.path.isdir(os.path.join(path, x))] #save foler/filename
    all_files.sort()

    info = pd.DataFrame()


    for f in all_files:
        
        if pattern == 'summaryinfo': 
            try: 
                df = pd.read_csv((path + '/' + f), 
                             header= 0,  sep ="|", engine = "python", names = ['ID', 'Job id', 'Hit', 'Confidence (%)',
                                                                            'Sequence identity (%)', 'Alignment coverage (%)', 'Hit info 1', 
                                                                            'Hit info 2', 'Hit info 3'])
    
                df['batch'] = f.split('/')[0]
                info = info.append(df, ignore_index =True)

        
            except Exception as e: print(e)
#            except: 
#                print('wrong column formating, try summaryinfo_extended')
            
        elif pattern == 'summaryinfo_extended':
            try:
                df = pd.read_csv((path + '/' + f), 
                   header=None, skiprows=0, comment= '#', error_bad_lines= False,  sep ="|", engine = "python", names = ['ID', 'Job id', 'Rank', 'Hit', 'Confidence (%)',
                                                                    'Sequence identity (%)', 'Qstart', 'Qend', 'Hit info 1',
                                                                     'Hit info 2', 'Hit info 3', 'Hit info 4', 'blank'])
                df.drop(df[df.Rank > 1].index, inplace=True)
                df.drop('blank', axis = 1, inplace=True)
                df['batch'] = f.split('/')[0]
                info = info.append(df, ignore_index =True)

            except Exception as e: print(e)
#            except:
#                print('wrong column formating, try summaryinfo')
            
        
        else:
            print('pattern error: use summaryinfo or summaryinfo_extended')
            
    return info

def import_sizes(path):
    #import model measurement files
    df = pd.read_csv(path, header= None, sep ="\t", names = ["Job id","x", "y", "z"])
    
    #clean Job id field
    df["Job id"] = df['Job id'].apply(lambda x: x.split('.')[0])
    df["Job id"] = df['Job id'].apply(lambda x: x.strip())
    
    return df
         
            
def clean_summaryinfo(df, pattern):
    
    #clean job id field
    df.loc[:,"Job id"] = df['Job id'].apply(lambda x: x.strip())
    
    #column order
    if pattern == 'summaryinfo': 
        col = ['ID link', 'Job id', 'Hit', 'Confidence (%)', 'Sequence identity (%)',
                       'Alignment coverage (%)', 'Hit info 1', 'Hit info 2', 'Hit info 3', 'ID', 
                       'start', 'end', 'ecd_len']
                
            
    elif pattern == 'summaryinfo_extended':
        col = ['ID link', 'Job id', 'Hit', 'Confidence (%)', 'Sequence identity (%)',
                       'Qstart', 'Qend', 'Hit info 1', 'Hit info 2', 'Hit info 3', 'ID', 
                       'start', 'end', 'ecd_len']
                
    else:
        print('pattern error: use summaryinfo or summaryinfo_extended')
                            
    
    #create and clean ID link field
    df.loc[:,'ID link'] = df['ID'].apply(lambda x: re.split('\[([^[\]]*)\]', str(x))[0])
    
    #clean ID field 
    df.loc[:,'ID'] = df['ID'].apply(lambda x: str(x).split('_')[0])
    
    #temp series to save ecd start:end info
    temp = df['ID'].apply(lambda x: re.split('\[([^[\]]*)\]', str(x))[1])
    
    #get start or end values
    df['start'] = temp.apply(lambda x: re.split(':', str(x))[0])
    df['end'] = temp.apply(lambda x: re.split(':', str(x))[1])
    
    #get ECD length
    df.loc[:,'ecd_len'] = df['end'].apply(lambda x: float(x)) - df['start'].apply(lambda x: float(x))
    
    #order/specifiy columns
    df = df[col]
    
    return df
    
    
def clean_summaryinfo_extended(df):
    
    #clean job id field
    df.loc[:,"Job id"] = df['Job id'].apply(lambda x: x.strip())
    
    #create and clean ID link field
    df.loc[:,'ID link'] = df['ID'].apply(lambda x: re.split('\[([^[\]]*)\]', str(x))[0])
    
    #clean ID field
    df.loc[:,'ID'] = df['ID'].apply(lambda x: str(x).split('_')[0])
    
    df.reset_index(drop=True, inplace = True) #need reset for skipped/droped lines
    
    #temp series to save ecd start:end info
    temp = df['ID'].apply(lambda x: re.split('\[([^[\]]*)\]', str(x))[1])
    
    #get start or end values
    df['start'] = temp.apply(lambda x: re.split(':', str(x))[0])
    df['end'] = temp.apply(lambda x: re.split(':', str(x))[1])
    
    #get ECD length
    df.loc[:,'ecd_len'] = df['end'].apply(lambda x: float(x)) - df['start'].apply(lambda x: float(x))
    
    #calculate Alignment coverage (%)
    df.loc[:,'Alignment coverage (%)'] = (df['Qend'] - df['Qstart']) / df['ecd_len']
    
    #column order
    col = ['ID link', 'Job id', 'Hit', 'Confidence (%)', 'Sequence identity (%)', 
           'Alignment coverage (%)', 'Qstart','Qend', 'ecd_len', 'Hit info 1', 
           'Hit info 2','Hit info 3', 'Hit info 4', 'ID', 'start', 'end']
    
    #order/specifiy columns
    df = df[col]
    
    return df
    
    
def find_dims(df):
    #find max and min dimensions, need dataframe with measured dims
    df['max_dim'] = df[["x", "y", "z"]].max(axis=1)
    df['min_dim'] = df[["x", "y", "z"]].min(axis=1)
    
    return df

def size_info_merge(info, sizes):
    #merge those dataframes!
    df = pd.merge(info, sizes, on = 'Job id', how ='outer')
    
    return df

def largest_unique_ecd(df):
    #select for the largest ecd when there are multiple measured ecds
    
    idx = df.groupby(['ID link'])['max_dim'].transform(max) == df['max_dim']
    
    return df[idx]

def largest_unique_ecd_seq(df):
    #select for the largest ecd when there are multiple measured ecds
    
    idx = df.groupby(['ID link'])['seq_len'].transform(max) == df['seq_len']
    
    return df[idx]

def well_modeled(df, confidence = 90, coverage = 70):
    #identify the models above a confidence of sequence converage threshold
    df = df[(df['Confidence (%)'] >= confidence) & (df['Alignment coverage (%)'] >= coverage)]
    
    return df
    

def surfacome_merge(info, surfaceome):
    #merge list all identified ECDs
    #need to used start and stop of ecd for the merge to specify specific ECDs
    info['start'] = info['start'].astype(float).round(0).astype(int)
    info['end'] = info['end'].astype(float).round(0).astype(int)
    
    surfaceome['start'] = surfaceome['start'].astype(float).round(0).astype(int)
    surfaceome['end'] = surfaceome['end'].astype(float).round(0).astype(int)
    
    df = pd.merge(info, surfaceome, on = ['ID link','start', 'end'], how ='outer')

    df = df.assign(ID = df['ID link']+'['+df['start'].astype(str)+':'+df['end'].astype(str)+']')
    
    return df 


#Testing code


#set path for protien measurements
#path = '/Users/southk/Box Sync/surfaceome_modeling/data/raw_data/models/UP000000437'
            
#file pattern
#pattern = 'summaryinfo'
#pattern = 'summaryinfo_extended'
#pattern = 'some random string'

#test!
#info = read_summaryinfo(path, pattern)
