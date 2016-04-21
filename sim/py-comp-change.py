import numpy as np
import MDAnalysis as mda
import sys

num = int(sys.argv[1]) #Must be a number. Number to be converted.
inf = sys.argv[2] #Must be a .gro file
outfile = sys.argv[3] #Must be a name, not a file with extension. Will be exported as gro.

# input_file = "dppc_bilayer.gro"
input_file = inf
bilayer = mda.Universe(input_file)

# Randomly select 35 lipids in range [1,129)
randomselection = np.random.randint(1,129,num)
rlist = randomselection.tolist()
sel = bilayer.select_atoms("resid 0") #This doesn't select anything, it just defines sel as an empty atom selection

for i in rlist:
    sel = sel + bilayer.select_atoms("resid %s" % i)

# Rename these lipids to DBPC
sel.resnames = "DBPC"

# Select all DPPC and W to transfer to output file
DPPW = bilayer.select_atoms("resname DPPC") + bilayer.select_atoms("resname W")

# Select the atoms from DBPC that actually belong to all DBPC molecules. IE, chop molecules by not selecting them.
rm = bilayer.select_atoms("resname DBPC and not (name C2A or name C3A or name C4A or name C2B or name C3B or name C4B)")

# Combine the unedited DPPC and W with the edited DBPC molecules. Then convert this to an atom list and write it to a gro file.
new = mda.Merge(rm, DPPW)
new = new.atoms
new.write(filename=outfile, format="GRO")