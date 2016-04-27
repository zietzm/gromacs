import sys
import MDAnalysis as mda

input_file = sys.argv[1] #GRO FILE
outfile = sys.argv[2] #FILE NAME NO EXTENSION

bilayer = mda.Universe(input_file)
DPPC = bilayer.select_atoms("resname DPPC")
DTPC = bilayer.select_atoms("resname DTPC")
W = bilayer.select_atoms("resname W")
selection = DPPC + DTPC + W
selection.write(filename=outfile, format="GRO")
