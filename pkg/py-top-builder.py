import sys
import MDAnalysis as mda
import numpy as np

infile = sys.argv[1] #GRO file
outtop = sys.argv[2] #FULL file name TOP
new = mda.Universe(infile)

# Count the number of each molecule.
numDBPC = new.select_atoms("resname DBPC")
numDPPC = new.select_atoms("resname DPPC")
numW = new.select_atoms("resname W")

#Account for the number of atoms in each molecule
ndb = len(np.unique(numDBPC.resids))
ndp = len(np.unique(numDPPC.resids))
nwat = len(np.unique(numW.resids))

#Comment out water if there isn't any
comm = ''
if nwat == 0:
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
