#!/bin/bash

### LIPID BICELLE ASSEMBLY
### BY MICHAEL ZIETZ
### COMPLETED 07 APR 2016

### REQUIRED FILES:
# dppc_bilayer.gro (FROM MARTINI WEB)
# minimization.mdp
# martini_md.mdp
# PDB_to_CSV.py
# Composition_Change.py
# TopologyBuilder.py
# Remove_Water.py
# ChangeBox.py
# water.gro
# dppc_single.itp
# dbpc_single.itp
# martini_v2.1.itp

clear
echo "Running Michael Zietz's Gromacs Simulation";

# cd vmBicelle6

#Define file names so python scripts can be applied generally.
InitGRO=dppc_bilayer.gro
InitPDB=dppc_bilayer.pdb
InitCSV=dppc_bilayer.csv
NewCompPDB=Out.pdb
NewCompCSV=Out.csv
NewCompTOP=topol.top
MixGRO=MixedBilayer.gro
MixMinTPR=MB_min.tpr
MixMin=MB_min
MixMinGRO=MB_min.gro
MixMinPDB=MB_min.pdb
MixMinCSV=MB_min.csv
MixNoWPDB=nowater_MB.pdb
MixNoWTOP=topol1.top
MixNoWGRO=nowater_MB.gro
MixNoWCSV=nowater_MB.csv
BigNoWGRO=nowaterReplicated.gro
BoxNoWGRO=BoxNoW.gro
SolvMBGRO=SolvMB.gro
SolvMBPDB=SolvMB.pdb
SolvMBCSV=SolvMB.csv
SolvMBTOP=topol2.top
SolvMinTPR=SolvMin.tpr

source /usr/local/gromacs/bin/GMXRC

#Convert Martini bilayer from gro to pdb to csv
editconf -f $InitGRO -o $InitPDB
python PDB_to_CSV.py $InitPDB $InitCSV

#Change csv composition. Revert output csv to pdb to gro. Output top
python Composition_Change.py $InitCSV $NewCompPDB $NewCompCSV
python TopologyBuilder.py $NewCompCSV $NewCompTOP
editconf -f $NewCompPDB -o $MixGRO

#Minimize the changed bilayer
grompp -f minimization.mdp -c $MixGRO -p $NewCompTOP -maxwarn 10 -o $MixMinTPR
mdrun -deffnm $MixMin -v -nt 1

#Remove water and output the new top gro pdb csv files
editconf -f $MixMinGRO -o $MixMinPDB
python PDB_to_CSV.py $MixMinPDB $MixMinCSV
python Remove_Water.py $MixMinCSV $MixNoWPDB
editconf -f $MixNoWPDB -o $MixNoWGRO

#Replicate the bilayer
genconf -f $MixNoWGRO -o $BigNoWGRO -nbox 3 3 1

#Change box size
python ChangeBox.py $BigNoWGRO $BoxNoWGRO

#Solvate
genbox -cp $BoxNoWGRO -cs water.gro -maxsol 20000 -vdwd 0.21 -o $SolvMBGRO

#Generate new topology
editconf -f $SolvMBGRO -o $SolvMBPDB
python PDB_to_CSV.py $SolvMBPDB $SolvMBCSV
python TopologyBuilder.py $SolvMBCSV $SolvMBTOP

#Minimize
grompp -f minimization.mdp -c $SolvMBGRO -p $SolvMBTOP -maxwarn 10 -o $SolvMinTPR
mdrun -deffnm SolvMin -v -nt 1

# Production Run
grompp -f martini_md.mdp -c SolvMin.gro -p $SolvMBTOP -maxwarn 10 -o SolvMB_Martini.tpr
tmux new-session -d -s Production_Martini_Run 'mdrun -deffnm SolvMB_Martini -v'
tmux detach -s Production_Martini_Run
