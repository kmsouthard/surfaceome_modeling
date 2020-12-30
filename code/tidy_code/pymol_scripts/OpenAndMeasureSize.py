#This script opens a pdb file and records the protein demensions as measured by the "Draw_Protein_Dimensions" script
#Draw_Protein_Dimensions was modified to print the dimensions to a file: currently listed as dist.txt
#may change to allow differing file names in the future, for now file name must be changed in "Draw_Protein_Dimensions.py"


from pymol import cmd

#clear pymol workspace
cmd.reinitialize(what="everything")

#select file, this will become redudented once file selection is automated
path = "/Users/southk/Box Sync/Fletcher Lab Data/Analysis/macrophage_surface_analysis/bmdm/"
#the unique part of the file name and part that is used within pymol
filename = "860c7a57fff6b17b"

#load file on path
cmd.load( path+filename+".final.pdb" )

#draw protein dimensions
cmd.run("/Users/southk/Box Sync/Fletcher Lab Data/Analysis/pymol/scripts/Draw_Protein_Dimensions.py")

#initialize file and append pdb filename (name must match file name in draw_Protein_Dimensions.py
f=open('dist.txt','a')
f.write(filename+"\t")
f.close()

# calculate the distance and store it in dst
draw_Protein_Dimensions(filename)

# write the formatted value of the distance (dst)
# to the output file
