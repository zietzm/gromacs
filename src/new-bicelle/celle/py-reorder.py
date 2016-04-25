import sys
import MDAnalysis as mda

input_file = sys.argv[1] #GRO FILE
outfile = sys.argv[2] #FILE NAME NO EXTENSION

bilayer = mda.Universe(input_file)
DXPC = bilayer.select_atoms("resname DXPC")
DTPC = bilayer.select_atoms("resname DTPC")
W = bilayer.select_atoms("resname W")
selection = DXPC + DTPC + W
selection.write(filename=outfile, format="GRO")
