import os
import fnmatch
import pandas as pd
import numpy as np

pattern = '*.pdb'
output_file = '200430_UP000005640_pdb_check.txt'

#select file, this will become redudented once file selection is automated
path = "/Users/southk/Box Sync/surfaceome_modeling/data/raw_data/models/UP000005640"
#the unique part of the file name and part that is used within pymol

#function for filename cleanup
def substract(a, b):                              
    return "".join(a.rsplit(b))

#calculate 3D distance between atoms
def dist_3d(row1, row2):
    #print(row1, row2)
    p1 = np.array([row1['x'], row1['y'], row1['z']])
    p2 = np.array([row2['x'], row2['y'], row2['z']])

    squared_dist = np.sum((p1-p2)**2, axis=0)
    try: 
        return np.sqrt(squared_dist).item()
    except:
        return np.nan


i = 0
for (path, dirs, files) in os.walk(path):
    for filename in fnmatch.filter(files, pattern):

        #print file name and fileID
        print(filename)
        fileID = substract(str(filename), ".pdb")
        print(fileID)

        #initialize file and append pdb filename (name must match file name in draw_Protein_Dimensions.py
        f=open(output_file,'a')
        f.write(fileID+"\t")
        f.close()
        
        #read in pdb 
        pdb_fw = pd.read_fwf(os.path.join(path, filename), header = None,
        skiprows = [0,1,2,3], widths =[5,6,5,4,6,12,8,8,6,6,12], 
        names= ['ATOM', 'ATOM_num', 'ATOM_id', 'resi_name', 'resi_num', 'x', 'y', 'z', 'one', 'zero', 'ATOM_name'], skipfooter = 1, engine = 'python')
        
        #remove whitespace
        pdb = pdb_fw.applymap(lambda x: x.strip() if type(x)==str else x)

        #limit dataframe to c-alphas
        pdb_ca = pdb[pdb['ATOM_id'] == 'CA'].reset_index()

        #calculate distance between C-alphas
        pdb_ca['dist'] = pdb_ca.apply(lambda x: dist_3d(pdb_ca[x.name:x.name+1:].reset_index(),pdb_ca[x.name+1:x.name+2:].reset_index()), axis = 1)
        ##initialize file and append pdb filename (name must match file name in draw_Protein_Dimensions.py

        #test for invalid distance
        test = pdb_ca[pdb_ca.dist > 5].empty
        max_dist = pdb_ca[pdb_ca.dist > 5]['dist'].max()
        count = pdb_ca[pdb_ca.dist > 5]['dist'].count()
        

        #write result to file
        f=open(output_file,'a')
        f.write(str(test)+"\t"+str(max_dist)+"\t"+str(count)+"\n")
        f.close()

                
        #                   
        #i += 1
        #
        #if i >= 10:
         #   break
