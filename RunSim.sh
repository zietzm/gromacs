#!/bin/bash
clear

echo "Running Michael Zietz's Gromacs Simulation";

cd vmBicelle6
source /usr/local/gromacs/bin/GMXRC

python Composition_Change.py

editconf -f Out.pdb -o MixedBilayer.gro

grompp -f minimization.mdp -c Out.gro -p topol.top -maxwarn 10 -o Out_min.tpr

mdrun -deffnm Out_min -v -nt 1

genconf -f Out_min.gro -o Out1.gro -nbox 3 3 1

# genbox -cp Out3.gro -cs water.gro -maxsol 18000 -vdwd 0.21 -o Out4.gro
#
# grompp -f minimization.mdp -c Out5.gro -p topol.top -maxwarn 10 -o Out5_min.tpr
#
# mdrun -deffnm Out5_min -v -nt 1
#
# grompp -f martini_md.mdp -c Out5_min.gro -p topol.top -maxwarn 10 -o Out5_Martini.tpr
#
# mdrun -deffnm Out5_Martini -v
