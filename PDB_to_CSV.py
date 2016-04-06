import numpy as np
import csv

# We want to convert a pdb file to a pure csv file.
infile = 'MB_min.pdb'

data = np.genfromtxt(infile, delimiter=" ", dtype=None)
