#openAndMeasureSizeMultiObject.py
#This script opens a pdb file and records the protein demensions as measured by the "Draw_Protein_Dimensions" script
#Draw_Protein_Dimensions was modified to print the dimensions to a file: currently listed as dist.txt
#may change to allow differing file names in the future, for now file name must be changed in "Draw_Protein_Dimensions.py"

import pymol
from pymol import cmd
import os
import fnmatch

pattern = '*.pdb'
output_file = '191019_dist_lig_receptor.txt'

#select file, this will become redudented once file selection is automated
path = "/Users/southk/Box Sync/Fletcher Lab Data/Analysis/macrophage_surface_analysis/bmdm/ligand receptor models/Ligand_Receptor_PDBs"
#the unique part of the file name and part that is used within pymol

#function for filename cleanup
def substract(a, b):                              
    return "".join(a.rsplit(b))

i = 0
for (path, dirs, files) in os.walk(path):
    for filename in fnmatch.filter(files, pattern):
        
        #clear pymol workspace
        cmd.reinitialize(what="everything")
        
        print(filename)
        fileID = substract(str(filename), ".pdb")
        print(fileID)
        
        #load file on path
        cmd.load(os.path.join(path, filename))
        
        #draw protein dimensions
        cmd.run("/Users/southk/Box Sync/surfaceome_modeling/code/raw_code/pymol_scripts/Draw_Protein_Dimensions.py")

        #initialize file and append pdb filename (name must match file name in draw_Protein_Dimensions.py
        f=open(output_file,'a')
        f.write(fileID+"\t")
        f.close()
        
        # calculate the distance and store it in dst
        try:
            draw_Protein_Dimensions('all', output_file)
        except KeyError:
            cmd.select('bad', 'name "hd*"')
            cmd.alter('bad', 'elem="H"')
            draw_Protein_Dimensions('all', output_file)
