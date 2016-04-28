import MDAnalysis as mda
import sys
infile=sys.argv[1]
outname=sys.argv[2]

bilayer = mda.Universe(infile)
nowaterselection = bilayer.select_atoms("not name W")
nowaterselection.write(filename=outname,format="GRO")