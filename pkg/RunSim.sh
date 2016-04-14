#!/bin/bash

### LIPID BICELLE ASSEMBLY
###
### LAST UPDATED 14 APR 2016

### REQUIRED FILES:
# dppc_bilayer.gro (FROM MARTINI WEB)
# minim.mdp
# martini_md.mdp
# water.gro
# dppc_single.itp
# dbpc_single.itp
# martini_v2.1.itp
# MDACompChange.py
# MDATopBuilder.py
# MDARemoveWater.py
# ChangeBox.py

source /usr/local/gromacs/bin/GMXRC

#Edit a number of lipids in the DPPC bilayer to be DBPC. Output top
python MDACompChange.py 35 dppc_bilayer.gro NewComp
python MDATopBuilder.py NewComp.gro NewComp

#Minimize the changed bilayer
grompp -f minim.mdp -c NewComp.gro -p NewComp.top -maxwarn 10 -o MinNC.tpr
mdrun -deffnm MinNC -v -nt 1

#Remove water and replicate bilayer
python MDARemoveWater.py MinNC.gro NoW
genconf -f NoW.gro -o BigNoW.gro -nbox 3 3 1

#Change box size
python ChangeBox.py BigNoW.gro BoxNoW.gro

#Solvate
genbox -cp BoxNoW.gro -cs water.gro -maxsol 31000 -vdwd 0.21 -o Solv.gro

#Generate new topology
python MDATopBuilder.py Solv.gro Solv

#Minimize
grompp -f minim.mdp -c Solv.gro -p Solv.top -maxwarn 10 -o SolvMin.tpr
mdrun -deffnm SolvMin -v -nt 1

#SYSTEM IS BLOWING UP BECAUSE ENERGIES TOO HIGH. NEED MORE PRELIM. MINIMIZATION
#INCLUDE POSITION RESTRAINT FILE IN TOPOLOGY FOR THIS MINIMIZATION.
# grompp -f MINIM.mdp -c $SolvMinGRO -p $SolvMBTOP -maxwarn 10 -o $SolvMin2TPR
# mdrun -deffnm $SolvMin2 -v

# Production Run
grompp -f martini_md.mdp -c SolvMin.gro -p Solv.top -maxwarn 10 -o SolvMartini.tpr
tmux new-session -d -s martini_run 'mdrun -deffnm SolvMartini -v'
tmux detach -s martini_run

# Sort all output files
mkdir inputs
mkdir results
mkdir intermediates
mkdir general-files
mv water.gro dppc_bilayer.gro dbpc_single.itp dppc_single.itp inputs/
mv NewComp.gro NewComp.top MinNC.tpr MinNC.log MinNC.gro MinNC.xtc NoW.gro BigNoW.gro BoxNoW.gro Solv.gro Solv.top SolvMin.tpr SolvMin.log SolvMin.xtc SolvMin.gro intermediates/
mv minim.mdp MDATopBuilder.py MDACompChange.py MDARemoveWater.py ChangeBox.py general-files/
mv SolvMartini.gro SolvMartini.log SolvMartini.tpr results/
