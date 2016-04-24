#!/bin/bash
source /usr/local/gromacs/bin/GMXRC

# Add 128 DXPC
genbox -ci str-dxpc-single.gro -box 6.31915 6.46099 10.05482 -nmol 128 -try 1000 -o gmx-vac-lipids.gro

#Solvate
genbox -cp gmx-vac-lipids.gro -cs str-water.gro -vdwd 0.21 -nmol 2000 -o gmx-lipids.gro

#Generate topology
python py-topology-builder.py gmx-lipids.gro gmx-lipids.top

#Energy minimize
grompp -f em.mdp -c gmx-lipids.gro -p gmx-lipids.top -maxwarn 10 -o em-lipids.gro

#Production run

