import numpy as np
import csv
import sys


infile = sys.argv[1]
outpdb = sys.argv[2]
outtop = sys.argv[3]

# infile = 'dppc_bilayer.csv'
# outpdb = 'Out.pdb'
# outtop = 'topol.top'


# We need 35% DBP. Since we already have 128 DPP, we will convert 44.8 ~ 45 or 35% to DBP.
number_to_be_converted = 45


data = np.genfromtxt(infile, delimiter=",", dtype=None)
ldat = data.size
counter0 = 0
dconc = number_to_be_converted*6 #We remove six beads to convert each DPPC to DBPC

#Remove the desired number of DPPC-unique beads to make DBPC.
for i in range(0,ldat):
    for j in ['C2A','C3A','C4A','C2B','C3B','C4B']:
        if counter0 < dconc:
            if (j in data[i]) == True:
                counter0 = counter0+1
                data = np.delete(data, (i))

#Change the edited rows of each molecule to read, DBP instead of DPP.
for k in range(0,dconc):
    if ('NC3' in data[k]) == True:
        for l in range(k,k+6):
            ((data[l])[3]) = 'DBP'

ldat = data.size # Update the size to reflect edits made.

# Print Output in PDB format. We could edit title, box size, etc.
ann=("TITLE     MIXED BILAYER" + '\n'
     "REMARK    THIS IS A SIMULATION BOX" +'\n'
     "CRYST1   63.191   64.610  100.548  90.00  90.00  90.00 P 1           1" +'\n'
     "MODEL        1" +'\n')
fmt="%0s%7s%5s%4s%6s%12s%8s%8s%6s%6s" # Print the tuples to be exactly spaced as pdb.
x=''
for i in range(0,ldat):
   x = x + (fmt % (tuple(data[i]))) + '\n' #format data as it should be for PDB

conc = ("TER"+'\n'
	"ENDMDL")

out = open(outpdb,'w') # Create an output file and print our annotations/data
out.write(ann)
out.write(x)
out.write(conc)
out.close()

#Generate topology file
counter1 = 0
counter2 = 0
counter3 = 0
for m in range(0,ldat):
    if ('DPP' in data[m]) == True:
        counter1 = counter1 + 1
    if ('DBP' in data[m]) == True:
        counter2 = counter2 + 1
    if ('W' in data[m]) == True:
        counter3 = counter3 +1
ndp = str(counter1/12) # Account for the number of beads in each molecule.
ndb = str(counter2/6)
nwat = str(counter3)

top = open(outtop, 'w')
top.write('#include "martini_v2.1.itp"'+ '\n'
           '#include "dppc_single.itp"'+'\n'
           '#include "dbpc_single.itp"'+'\n'
           '\n'
          '[ system ]' +'\n'
          'LIPID BICELLE SELF-ASSEMBLY'+'\n'
          '\n'
          '[ molecules ]'+'\n'
          'DPPC %s' % ndp +'\n'
          'DBPC %s' % ndb +'\n'
          'W %s' % nwat +'\n')
top.close()
