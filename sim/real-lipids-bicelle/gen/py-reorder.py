import sys
import MDAnalysis as mda

input_file = sys.argv[1] #GRO FILE
outfile = sys.argv[2] #FILE NAME NO EXTENSION

bilayer = mda.Universe(input_file)
DBPC = bilayer.select_atoms("resname DBPC")
DPPC = bilayer.select_atoms("resname DPPC")
W = bilayer.select_atoms("resname W")
selection = DBPC + DPPC + W
selection.write(filename=outfile, format="GRO")
