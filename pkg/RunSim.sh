
#!/bin/bash

### LIPID BICELLE ASSEMBLY
###
### LAST UPDATED 09 APR 2016

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
echo "Running Gromacs Simulation";

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
SolvMin=SolvMin
SolvMinGRO=SolvMin.gro
SolvMinTPR=SolvMin.tpr
SolvMin2TPR=SolvMin2.tpr
SolvMin2=SolvMin2
SolvMin2GRO=SolvMin2.gro
SolvMartiniTPR=SolvMartini.tpr
SolvMartini=SolvMartini3


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
python Remove_Water.py $MixMinCSV $MixNoWPDB #### ERROR HERE IN PY
editconf -f $MixNoWPDB -o $MixNoWGRO

#Replicate the bilayer
genconf -f $MixNoWGRO -o $BigNoWGRO -nbox 3 3 1

#Change box size
python ChangeBox.py $BigNoWGRO $BoxNoWGRO

#Solvate
genbox -cp $BoxNoWGRO -cs water.gro -maxsol 31000 -vdwd 0.21 -o $SolvMBGRO

#Generate new topology
editconf -f $SolvMBGRO -o $SolvMBPDB
python PDB_to_CSV.py $SolvMBPDB $SolvMBCSV
python TopologyBuilder.py $SolvMBCSV $SolvMBTOP

#Minimize
grompp -f minimization.mdp -c $SolvMBGRO -p $SolvMBTOP -maxwarn 10 -o $SolvMinTPR
mdrun -deffnm $SolvMin -v -nt 1

#SYSTEM IS BLOWING UP BECAUSE ENERGIES TOO HIGH. NEED MORE PRELIM. MINIMIZATION
# grompp -f MINIM.mdp -c $SolvMinGRO -p $SolvMBTOP -maxwarn 10 -o $SolvMin2TPR
# mdrun -deffnm $SolvMin2 -v

# Production Run
grompp -f martini_md.mdp -c $SolvMinGRO -p $SolvMBTOP -maxwarn 10 -o $SolvMartiniTPR
tmux new-session -d -s martini_run 'mdrun -deffnm SolvMartini -v'
tmux detach -s martini_run

# Sort all output files
mkdir inputs
mkdir results
mkdir intermediates
mkdir general-files
mv $InitGRO $InitPDB $InitCSV dbpc_single.itp dppc_single.itp inputs/
mv $NewCompPDB $NewCompCSV $NewCompTOP $MixGRO $MixMinTPR $MixMinGRO $MixMinPDB $MixMinCSV $MixNoWPDB $MixNoWTOP $MixNoWGRO $MixNoWCSV $BigNoWGRO $BoxNoWGRO $SolvMBGRO $SolvMBPDB $SolvMBCSV $SolvMinGRO $SolvMinTPR $SolvMin2TPR $SolvMin2GRO $SolvMartiniTPR intermediates/
mv cg_bonds.tcl ChangeBox.py Composition_Change.py Remove_Water.py PDB_to_CSV.py TopologyBuilder.py water.gro general-files/
mv $SolvMBTOP results/
