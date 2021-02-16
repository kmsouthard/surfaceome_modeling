#!/usr/bin/env python2
'''
Mobidb-lite disorder annotations
file could be processed into gff like file? and consensus disorder regions

takes proteome ID and finds disordered protein consensus annotations

Proteome ID
*UP000000803
*UP000005640
*UP000000589
*UP000002195
*UP000001940
*UP000000437

* = done


'''



import requests
import pandas as pd

# Define request

proteome = 'UP000000437'

surfaceome_path = '/Users/southk/Box Sync/surfaceome_modeling/data/raw_data/surface_proteins/'
output_path = '/Users/southk/Box Sync/surfaceome_modeling/data/raw_data/disorder_annotations/IUPred-long'

surfaceome = pd.read_csv(surfaceome_path+proteome+'.csv')

accessions = surfaceome['ID link'].unique().tolist()

acceptHeader = 'text/plain'
f= open(output_path+proteome+'_disorder.txt', 'w')

for accession in accessions:

    #create individual file name for proteome

    print("retrieving "+accession+" data...")


    url = "http://mobidb.bio.unipd.it/ws/"+accession+"/consensus"

    url = requests.get(url, headers={"Accept" : acceptHeader})

    data = url.text

    for line in data.splitlines():
        print(line)
        if 'mobidb_consensus.disorder.predictors.mobidb-lite.regions' in line:
            #print(line)
            f.write(line+'\n')

f.close()


# handle data
#print (data)
