import MDAnalysis as mda
import sys
import random

input_file = sys.argv[1] #GRO FILE
outfile = sys.argv[2] #FILE NAME NO EXTENSION

bilayer = mda.Universe(input_file)
b = bilayer.atoms
selection = b.select_atoms("resid 0")
rlist = random.sample(xrange(1, 129), 35)
for i in rlist:
    selection = selection + bilayer.select_atoms("resid %s" % i)
selection.resnames = "DBPC"
selection = bilayer.select_atoms("resname DBPC and not (name C2A or name C3A or name C4A or name C2B or name C3B or name C4B)")
selection = selection + bilayer.select_atoms("resname DPPC or resname W")
selection.write(filename=outfile, format="GRO")
