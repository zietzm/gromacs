#!/bin/bash
clear

echo "Running Michael Zietz's Gromacs Simulation";

# cd vmBicelle6
source /usr/local/gromacs/bin/GMXRC

FILE0=dppc_bilayer.gro
FILE1=dppc_bilayer.pdb
FILE2=dppc_bilayer.csv
FILE3=Out.pdb
FILE4=topol.top
FILE5=MB_min.gro
FILE6=MB_min.pdb
FILE7=MB_min.csv
FILE8=nowater_MB_min.pdb
FILE9=topol1.top
FILE10=Out1.gro
FILE11=Out2.gro
FILE12=Out3.pdb
FILE13=Out3.csv
FILE14=topol2.top
FILE15=



#Convert Martini bilayer from gro to pdb to csv
editconf -f dppc_bilayer.gro -o dppc_bilayer.pdb
python PDB_to_CSV.py $FILE1 $FILE2

#Change composition, change csv to pdb to gro and output top
python Composition_Change.py $FILE2 $FILE3 $FILE4
editconf -f Out.pdb -o MixedBilayer.gro

#PRELIMINARY MINIMIZATION
grompp -f minimization.mdp -c MixedBilayer.gro -p topol.top -maxwarn 10 -o MB_min.tpr
mdrun -deffnm MB_min -v -nt 1

#REMOVE WATER. Output pdb to gro and top
editconf -f MB_min.gro -o MB_min.pdb
python PDB_to_CSV.py $FILE6 $FILE7
python Remove_Water.py $FILE7 $FILE8 $FILE9
editconf -f nowater_MB_min.pdb -o nowater_MB_min.gro

#Replicate the bilayer
genconf -f nowater_MB_min.gro -o Out1.gro -nbox 3 3 1

#Change box size [MAKE A NEW PYTHON FILE]
python ChangeBox.py $FILE10 $FILE11

#Solvate
genbox -cp Out2.gro -cs water.gro -maxsol 20000 -vdwd 0.21 -o Out3.gro

#Generate new topology
editconf -f Out3.gro -o Out3.pdb
python PDB_to_CSV.py $FILE12 $FILE13
python TopologyBuilder.py $FILE13 $FILE14

#Minimize
grompp -f minimization.mdp -c Out3.gro -p topol2.top -maxwarn 10 -o Out3_min.gro
mdrun -deffnm Out3_min -v -nt 1

# Production Run
grompp -f martini_md.mdp -c Out3_min.gro -p topol2.top -maxwarn 10 -o Out3_Martini.tpr
tmux
mdrun -deffnm Out3_Martini -v
