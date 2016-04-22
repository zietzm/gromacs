We have started using functional prefixes for each file to make things 
more simple and clear.

The prefixes are as follows:
* em- energy minimization files
* gmx- general gromacs output files from python, etc. not mdrun.
* md- bigger molecular dynamics run files
* py- python files
* str- initial input structure files
* top- topology input files. Both structure and position restraints

Required files to begin:
cg_bonds.tcl
em-max.mdp
md-martini.mdp
md-posre.mdp
py-change-box.py
py-comp-change.py
py-posre-all-top-builder.py
py-posre-top-builder.py
py-remove-water.py
py-top-builder.py
README.md
Run.sh
str-dppc-bilayer.gro
str-water.gro
top-dbpc-posre.itp
top-dbpc-posre-all.itp
top-dbpc-single.itp
top-dppc-posre.itp
top-dppc-posre-all.itp
top-dppc-single.itp
top-martini-v2.1.itp
