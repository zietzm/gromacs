#!/bin/bash
source /usr/local/gromacs/bin/GMXRC

# Add 128 DXPC
genbox -ci str-dxpc-single.gro -box 6.31915 6.46099 10.05482 -nmol 128 -try 1000 -o gmx-vac-lipids.gro

#Generate topology
python py-top-builder.py gmx-vac-lipids.gro gmx-vac.top

#Energy minimize
grompp -f em.mdp -c gmx-vac-lipids.gro -p gmx-vac-top.top -maxwarn 10 -o em-vac.tpr
mdrun -deffnm em-vac -v

#Solvate and increase box size to be enough for 2000 lipids
genbox -cp em-vac.gro -cs str-water.gro -vdwd 0.21 -maxsol 2000 -box 6.5 6.5 12.5 -o gmx-lipids.gro

#Generate topology
python py-top-builder.py gmx-lipids.gro gmx-lipids.top

#Minimization
grompp -f em.mdp -c gmx-lipids.gro -p gmx-lipids.top -maxwarn 10 -o em-lipids.tpr
mdrun -deffnm em-lipids -v

#Production run
grompp -f md-bilayer.mdp -c em-lipids.gro -p gmx-lipids.top -maxwarn 10 -o md-lipids.tpr
mdrun -deffnm md-lipids -v

#Put all molecules into the frame
trjconv -f md-lipids.xtc -s md-lipids.tpr -pbc mol -o gmx-bilayer.xtc
trjconv -f gmx-bilayer.xtc -s md-lipids.tpr -dump 30000 -o gmx-bilayer.gro

#Delete Unnecessary files and reorganize
mv gmx-lipids.top gmx-bilayer.top
rm em-lipids.edr em-lipids.gro em-lipids.log em-lipids.tpr tm-lipids.trr gmx-lipids.gro gmx-lipids.xtc md-lipids.cpt md-lipids.edr md-lipids.gro md-lipids.log md-lipids.tpr md-lipids.xtc mdout.mdp 
