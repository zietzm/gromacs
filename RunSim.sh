#!/bin/bash
clear

echo "Running Michael Zietz's Gromacs Simulation";

# cd vmBicelle6
source /usr/local/gromacs/bin/GMXRC

#Define file names so python scripts can be applied generally.
InitGRO=dppc_bilayer.gro
InitPDB=dppc_bilayer.pdb
InitCSV=dppc_bilayer.csv
NewCompPDB=Out.pdb
NewCompTOP=topol.top
MixMinGRO=MB_min.gro
MixMinPDB=MB_min.pdb
MixMinCSV=MB_min.csv
MixNoWPDB=nowater_MB_min.pdb
MixNoWTOP=topol1.top
BigNoWGRO=nowaterReplicated.gro
BoxNoWGRO=BoxNoW.gro
SolvMBPDB=SolvMB.pdb
SolvMBCSV=SolvMB.csv
SolvMBTOP=topol2.top



#Convert Martini bilayer from gro to pdb to csv
editconf -f dppc_bilayer.gro -o dppc_bilayer.pdb
python PDB_to_CSV.py $InitPDB $InitCSV

#Change composition, change csv to pdb to gro and output top
python Composition_Change.py $InitCSV $NewCompPDB $NewCompTOP
editconf -f Out.pdb -o MixedBilayer.gro

#PRELIMINARY MINIMIZATION
grompp -f minimization.mdp -c MixedBilayer.gro -p topol.top -maxwarn 10 -o MB_min.tpr
mdrun -deffnm MB_min -v -nt 1

#REMOVE WATER. Output pdb to gro and top
editconf -f MB_min.gro -o MB_min.pdb
python PDB_to_CSV.py $MixMinPDB $MixMinCSV
python Remove_Water.py $MixMinCSV $MixNoWPDB $MixNoWTOP
editconf -f nowater_MB_min.pdb -o nowater_MB_min.gro

#Replicate the bilayer
genconf -f nowater_MB_min.gro -o nowaterReplicated.gro -nbox 3 3 1

#Change box size
python ChangeBox.py $BigNoWGRO $BoxNoWGRO

#Solvate
genbox -cp BoxNoW.gro -cs water.gro -maxsol 20000 -vdwd 0.21 -o SolvMB.gro

#Generate new topology
editconf -f SolvMB.gro -o SolvMB.pdb
python PDB_to_CSV.py $SolvMBPDB $SolvMBCSV
python TopologyBuilder.py $SolvMBCSV $SolvMBTOP

#Minimize
grompp -f minimization.mdp -c SolvMB.gro -p topol2.top -maxwarn 10 -o SolvMB_min.gro
mdrun -deffnm SolvMB_min -v -nt 1

# Production Run
grompp -f martini_md.mdp -c SolvMB_min.gro -p topol2.top -maxwarn 10 -o SolvMB_Martini.tpr
tmux
mdrun -deffnm SolvMB_Martini -v
