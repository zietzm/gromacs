import MDAnalysis as mda
import sys
import random

input_file = sys.argv[1] #Must be a .gro file
outfile = sys.argv[2] #Must be a name, not a file with extension. Will be exported as gro.
# input_file = "dppc_bilayer.gro"

bilayer = mda.Universe(input_file)

# Randomly select 35 lipids in range [1,129)
rlist = rlist = random.sample(xrange(1, 129), 35)

sel = bilayer.select_atoms("resid 0") #This doesn't select anything, it just defines sel as an empty atom selection
for i in rlist:
    sel = sel + bilayer.select_atoms("resid %s" % i)

# Rename these lipids to DBPC
sel.resnames = "DBPC"

# Select all DPPC and W to transfer to output file
DPPW = bilayer.select_atoms("resname DPPC or resname W")

# Select the atoms from DBPC that actually belong to all DBPC molecules. IE, chop molecules by not selecting them.
rm = bilayer.select_atoms("resname DBPC and not (name C2A or name C3A or name C4A or name C2B or name C3B or name C4B)")

# Combine the unedited DPPC and W with the edited DBPC molecules. Then convert this to an atom list and write it to a gro file.
new = mda.Merge(rm, DPPW)
new = new.atoms
new.write(filename=outfile, format="GRO")

# Change the box dimensions automatically to be the same as the original bilayer.
file1 = open(str(outfile + '.gro'))
data=file1.read()

# Get dimensions from input file. We don't want to change those yet
df = open(input_file).read()
dimens = df[-31:]

data = data.replace("   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000",
dimens)

out = open(str(outfile + '.gro'),'w')
out.write(data)
out.close()
