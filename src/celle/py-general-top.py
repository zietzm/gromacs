# Note: usage in terminal is '$ python py-general-top.py input.gro output.top type'
# Where 'type' could be 'single', 'posre', or 'posre-all'. Posre is just lipid headgroups.

import sys
import MDAnalysis as mda
import numpy as np

infile = sys.argv[1]
outtop = sys.argv[2]
type_of_top = sys.argv[3] #either "single", "posre", or "posre-all"
new = mda.Universe(infile)

types = ["single","posre","posre-all"]

if (type_of_top in types) != True:
    raise ValueError("You must give a valid type of topology file")

resnames = list(np.unique(list(new.atoms.resnames)))
molecule_out = ''
itp_include = ''
for i in resnames:
    molecule_out = molecule_out +str(i + ' ' + 
        str(len(np.unique(new.select_atoms("resname %s" % i).resids))) + '\n')
    if i != 'W':
        itp_include = itp_include + str('#include "top-' + i.lower() + '-' + 
            type_of_top + '.itp"' + '\n')

top = open(outtop, 'w')
top.write(';::\n;;; Generated using py-general-top.py \n;;;\n\n'
          '#include "top-martini-v2.1.itp"'+ '\n'
          +itp_include+
          '\n'
          '[ system ]' +'\n'
          'LIPID SYSTEM TOPOLOGY'+'\n'
          '\n'
          '[ molecules ]'+'\n'
          + molecule_out)
top.close()
