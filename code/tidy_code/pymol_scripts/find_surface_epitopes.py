#This script opens a list of epitopes, finds the sequence, and displays the surface
#structure: 6CRZ

import pymol
from pymol import cmd
import os
import fnmatch


epitopes_file = open('/Users/southk/Box Sync/coronavirus_struct_exploration/epitopes/epitopes.txt', 'r')

output_file = 'epitopes_len.txt'

i = 0

for line in epitopes_file.readlines()[1:]:
    peptide = line
    print peptide

    i += 1
    findseq(needle = peptide, haystack= "6CRZ", selName= "sele_"+peptide, firstOnly=1)
    #findseq(peptide,'6CRZ', selName = peptide)
    listselection(peptide+'and chain A')
    print counter
    #f=open(output_file,'a')
    #f.write(values[3]+'\t'+values[1]+"\t"+values[2]+'\t'+values[4]+'\t')
    #f.close()


epitopes_file.close()

##function for filename cleanup
#def substract(a, b):                              
#    return "".join(a.rsplit(b))
#
#i = 0
#for index, row in domains.iterrows():