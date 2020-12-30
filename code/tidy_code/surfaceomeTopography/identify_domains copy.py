#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:25:46 2019

@author: southk

extracellular domain identification
    steps: load proteome gff and search for ECD annotation in 'notes' column
    how to deal with gpi anchored protiens: gpi/lypidation identified in gff
        some cases maybe not properly annotated eg CEACAM5 in mouse
        
for mouse proteome more annotations for gpi (150 vs 126) if tab file and
subcellular location [CC] column (GO term annotation) is used for ID
may need retrive annotations from both tab and gff files for full coverage
as proteomes may not all be annotated in the same way
    
    
unannotated proteins (uniprot) maybe need another identification method
cell membrane go terms and transmembrane domain annotation? Calculate ECD
from TMD annotation?
    
    
starting with just proteins that are annotated with an extracellular domain
    and gpi
    then compare to published surfaceome datasets, what is missing/how are
    they annotated
    Results: generally good human surfaceome coverage, most agreement with the
    machine learning surfaceome (2542 common/2886), both missed proteins and
    unique proteins their localization is uncertain

"""
import pandas as pd


def seq_len(df_gff):
    return df_gff.assign(seq_len = df_gff['end']-df_gff['start'])


def subset_ecds(df_gff):
    
    #returns subset of sequences annotated as Extracellular in gff file
    df = df_gff[df_gff['note'].str.contains('Extracellular', na=False)].reset_index(drop=True)
    return seq_len(df)

def subset_icds(df_gff):
    
    #returns subset of sequences annotated as cytoplasmic in gff file
    df = df_gff[df_gff['note'].str.contains('Cytoplasmic', na=False)].reset_index(drop=True)
    return seq_len(df)

def subset_e_and_i(df_gff):

	#returns annotated ecds and icds to look for unexposed n-termini
    df = df_gff[df_gff['note'].str.contains('Cytoplasmic|Extracellular', na=False)].reset_index(drop=True)
    return seq_len(df)	

def subset_gpi(df_tab, df_gff):
    
    #returns subset of sequences with cell membrane and gpi annotations (tab)
    #retrieves sequences for modeling from gff file
    gpi = df_tab.loc[
                df_tab['Subcellular location [CC]'].str.contains('Cell membrane|Membrane|cell membrane|membrane',na=False) 
              & 
                df_tab['Subcellular location [CC]'].str.contains('GPI',na=False)
                ]
    gpi_prots = df_gff.merge(gpi['Entry'], right_on = 'Entry', left_on = 'ID link', how = 'right')
    
    df = gpi_prots.loc[(gpi_prots['td'] == 'Chain') | (gpi_prots['td'] == 'Peptide')]
    return seq_len(df)  

def split_seq(seq_len):
    #split sequence in the middle
    first_end = round(seq_len / 2)
    second_start = first_end + 1    
    return first_end, second_start

def shorten_end(seq_len, start, end, c):
    first_end, second_start = split_seq(seq_len)
    if c == 0:
        return first_end + start
    if c == 1:
        return end
    
def truncate_start(seq_len, start, end, c):
    first_end, second_start = split_seq(seq_len)
    if c == 0:
        return start
    if c == 1:
        return (second_start+start)
    
def split_long(df):
    long = df[df['seq_len'] >= 3000]
    dup = pd.concat([long]*2, ignore_index=True)
    c = dup.groupby(['ID link']).cumcount()
    grouped = dup.assign(g = c)
    grouped['end'] = grouped.apply(lambda row: shorten_end(row['seq_len'], row['start'], row['end'], row['g']), axis=1)
    grouped['start'] = grouped.apply(lambda row: truncate_start(row['seq_len'], row['start'], row['end'], row['g']), axis=1)
    return grouped.dropna(subset=['start','end'])

def phyre_filter(df):
    phyre = df[(df['seq_len'] >= 30) & (df['seq_len'] < 3000)]
    return phyre.reset_index(drop=True)


def surfaceome_df(gpi, ecd, short, long = pd.DataFrame()):
    #combine dataframes for creating total ECDs dataframe
    df1 = gpi.assign(source = 'gpi')
    df2 = ecd.assign(source = 'ecd')
    df4 =short.assign(source = 'short')
    if long.empty == True:
        total = pd.concat([df1, df2, df4], sort = False)
        total.drop(['index', 'Entry'], axis =1, inplace = True)
        
    else:
        df3 =long.assign(source = 'long')
        total = pd.concat([df1, df2, df3, df4], sort = False)
        total.drop(['index', 'Entry', 'g'], axis =1, inplace = True)
        
    return total
    
