#This script fetchs a pdb file, selects the chain and residues, and records the domain demensions as measured by the "Draw_Protein_Dimensions" script
#Draw_Protein_Dimensions was modified to print the dimensions to a file: currently listed as dist.txt
#may change to allow differing file names in the future, for now file name must be changed in "Draw_Protein_Dimensions.py"

import pymol
from pymol import cmd
import os
import fnmatch

domains_file = open('/Users/southk/Box Sync/surfaceome_modeling/code/raw_code/human_surfaceome_pdbs.csv', 'r')

output_file = 'structure_sizes_200710.txt'

i = 0

for line in domains_file.readlines()[1:]:
    values = line.split(',')
    print values[1]

    i += 1
    #clear pymol workspace
    cmd.reinitialize(what="everything") 
    #load file on path
    cmd.fetch(values[1]+values[2])
    print values[1]
    print values[2]
    cmd.remove('solvent')
    #print 'resi ' +values[4]+ ' in chain '+values[2]
    #print 'resi %s in chain %s' %(values[4],values[2])

    chain_length = len(cmd.get_model(values[1]).get_residues())
    print len(cmd.get_model(values[1]).get_residues())

    cmd.select('dom', str(values[1]))
        
    #draw protein dimensions
    cmd.run("/Users/southk/Box Sync/surfaceome_modeling/code/raw_code/pymol_scripts/Draw_Protein_Dimensions.py")

        #initialize file and append pdb filename (name must match file name in draw_Protein_Dimensions.py
    f=open(output_file,'a')
    f.write(str.strip(values[2])+'\t'+values[1]+'\t'+str(chain_length)+'\t')
    f.close()
        
    # calculate the distance and store it in dst
    try:
        draw_Protein_Dimensions('dom', output_file)

    except KeyError:
        cmd.select('bad', 'name "hd*"')
        cmd.alter('bad', 'elem="H"')
        cmd.select('bad', 'name "d*"')
        cmd.alter('bad', 'elem="H"')
        draw_Protein_Dimensions('dom', output_file)
        
    except ZeroDivisionError:
        f=open(output_file,'a')
        f.write('\t\n')
        f.close()
                
    #if i > 10:
     #   break

domains_file.close()

##function for filename cleanup
#def substract(a, b):                              
#    return "".join(a.rsplit(b))
#
#i = 0
#for index, row in domains.iterrows():