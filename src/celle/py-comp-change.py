import MDAnalysis as mda
import sys
import random

input_file = sys.argv[1] #GRO FILE
outfile = sys.argv[2] #FILE NAME NO EXTENSION

bilayer = mda.Universe(input_file)
selection = bilayer.select_atoms("resid 0")

rlist = random.sample(xrange(1, 129), 35)
for i in rlist:
    selection = selection + bilayer.select_atoms("resid %s" % i)
    
selection.resnames = "DTPC"
selection = bilayer.select_atoms("resname DTPC and not (name C3A or name C4A or name C5A or name C6A or name C3B or name C4B or name C5B or name C6B)")
selection = selection + bilayer.select_atoms("resname DXPC or resname W")

selection.write(filename=outfile, format="GRO")
