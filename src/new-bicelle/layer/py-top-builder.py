import sys
import MDAnalysis as mda
import numpy as np

infile = sys.argv[1] #GRO file
outtop = sys.argv[2] #FULL file name TOP
new = mda.Universe(infile)

# Count the number of each molecule.
numDXPC = new.select_atoms("resname DXPC")
numDTPC = new.select_atoms("resname DTPC")
numW = new.select_atoms("resname W")

#Account for the number of atoms in each molecule
ndx = len(np.unique(numDXPC.resids))
ndt = len(np.unique(numDTPC.resids))
nwat = len(np.unique(numW.resids))

#Comment out water if there isn't any
comX = ''
if ndx == 0:
	comX = comX + ';'
comT = ''
if ndt == 0:
	comT = comT + ';'
comW = ''
if nwat == 0:
	comW = comW + ';'

top = open(outtop, 'w')
top.write('#include "top-martini-v2.1.itp"'+ '\n'
           '#include "top-dxpc-single.itp"'+'\n'
           '#include "top-dtpc-single.itp"'+'\n'
           '\n'
          '[ system ]' +'\n'
          'MIXED BILAYER'+'\n'
          '\n'
          '[ molecules ]'+'\n'
          '%sDXPC %s' % (comX, ndx) +'\n'
          '%sDTPC %s' % (comT, ndt) +'\n'
          '%sW %s' % (comW, nwat) +'\n')
top.close()
