import numpy as np
import sys

infile = sys.argv[1] #must be csv
outtop = sys.argv[2]

# infile = 'dppc_bilayer.csv'
# outpdb = 'Out.pdb'
# outtop = 'topol.top'


data = np.genfromtxt(infile, delimiter=",", dtype=None)

ldat = data.size

#Count numbers of each molecule
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

#Comment out water if there isn't any.
comm = ''
if counter3 == 0:
    comm = comm + ';'

top = open(outtop, 'w')
top.write('#include "martini_v2.1.itp"'+ '\n'
           '#include "dppc_single.itp"'+'\n'
           '#include "dbpc_single.itp"'+'\n'
           '\n'
          '[ system ]' +'\n'
          'MIXED BILAYER'+'\n'
          '\n'
          '[ molecules ]'+'\n'
          'DPPC %s' % ndp +'\n'
          'DBPC %s' % ndb +'\n'
          '%sW %s' % (comm, nwat) +'\n')
top.close()
