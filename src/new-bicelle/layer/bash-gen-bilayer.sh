#!/bin/bash
set -e

### Requires the following files:
### bash-gen-bilayer.sh
### em.mdp
### md-bilayer.mdp
### py-top-builder.py
### README.md
### str-dxpc-single.gro
### str-water.gro
### top-dtpc-single.itp
### top-dxpc-single.itp
### top-martini-v2.1.itp

source /usr/local/gromacs/bin/GMXRC
genbox -ci str-dxpc-single.gro -box 6.31915 6.46099 10.05482 -nmol 128 -try 1000 -o gmx-vac-lipids.gro
python py-top-builder.py gmx-vac-lipids.gro gmx-vac.top
grompp -f em.mdp -c gmx-vac-lipids.gro -p gmx-vac.top -maxwarn 10 -o em-vac.tpr
mdrun -deffnm em-vac -v
genbox -cp em-vac.gro -cs str-water.gro -vdwd 0.21 -maxsol 2000 -box 6.5 6.5 12.5 -o gmx-lipids.gro
python py-top-builder.py gmx-lipids.gro gmx-lipids.top
grompp -f em.mdp -c gmx-lipids.gro -p gmx-lipids.top -maxwarn 10 -o em-lipids.tpr
mdrun -deffnm em-lipids -v
grompp -f md-bilayer.mdp -c em-lipids.gro -p gmx-lipids.top -maxwarn 10 -o md-lipids.tpr
mdrun -deffnm md-lipids -v
echo "0" | trjconv -f md-lipids.xtc -s md-lipids.tpr -pbc mol -o gmx-bilayer.xtc
echo "0" | trjconv -f gmx-bilayer.xtc -s md-lipids.tpr -dump 30000 -vel -o gmx-bilayer.gro
cp gmx-lipids.top gmx-bilayer.top
