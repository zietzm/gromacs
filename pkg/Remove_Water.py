import numpy as np
import sys

infile = sys.argv[1] #needs to be a csv file
outpdb = sys.argv[2]

# infile = 'MB_min.csv'
# outpdb = 'nowater_MB_min.pdb'

data = np.genfromtxt(infile, delimiter=",", dtype=None)
ldat = data.size
   
#Count number of water
counter = 0
for m in range(0,ldat):
    if ('W' in data[m]) == True:
        counter = counter+1

num = int(data.size) - counter

#Remove that number of water
while (int(data.size) > 1536):
    data = np.delete(data,int(data.size - 1),0)

ldat = data.size # Update the number of atoms to reflect edits made.

# Print Output in PDB format. We could edit title, box size, etc.
ann=("TITLE     WATERLESS MIXED BILAYER" + '\n'
     "REMARK    THIS IS A SIMULATION BOX" +'\n'
     "CRYST1   63.191   64.610  100.548  90.00  90.00  90.00 P 1           1" +'\n'
     "MODEL        1" +'\n')
fmt="%0s%7s%5s%4s%6s%12s%8s%8s%6s%6s" # Print the tuples to be exactly spaced as pdb.
x=''
for i in range(0,ldat):
   x = x + (fmt % (tuple(data[i]))) + '\n' #We should print automatically to PDB

conc=("TER" + '\n'
      "ENDMDL")

out = open(outpdb,'w') # Create an output file and print our annotations/data
out.write(ann)
out.write(x)
out.write(conc)
out.close()
