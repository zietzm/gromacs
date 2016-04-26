# "Real" Lipid simulation

In the dissertation whose method I initially replicated to generate a bicelle, Jiang et al discuss 
DPPC and DBPC lipids. Their given definition of DPPC matches the DPPC given on the Martini site.
However, they describe DBPC as having one cg bead in each tail. Martini's definition of "DBPC" is 
a lipid with five cg beads in each tail. Because of this discrepancy, I will make a new simulation 
using DTPC and DXPC-- lipids of 2 and 6 cg tail beads, respectively.

This will require a large overhaul of many python files and input parameters.

We also have to create and minimize our own DXPC bilayer. In the previous run, we were able to use a
pre-made bilayer from Martini.

Outline:
1. Generate a bilayer using the 'layer' directory
2. Take the output of 'layer' and, using the 'celle' directory, convert it to a bicelle.

Specifics for each step can be found in subdirectory README files. It is the goal of this file structure
for bash files to be all that is required for total completion of each step.