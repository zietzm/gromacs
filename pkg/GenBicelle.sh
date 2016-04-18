#!/bin/bash

#Change composition
python MDAChangeComp.py dppc_bilayer.gro py-mixed-bilayer

#Output topology
python MDATopBuilder.py py-mixed-bilayer.gro py-mixed-top.top

#Energy minimize
grompp -f em-min.mdp -c py-mixed-bilayer.gro -p mixed-top.top -maxwarn 10 -o em-mix-bil.tpr
mdrun -deffnm em-mix-bil -v ##can we stop using deffnm and output only a few files

###___Remove water___###

#Replicate bilayer
genconf -f em-mix-bil.gro -o gmx-replicated.gro -nbox 3 3 1

#Change box size
python ChangeBox.py gmx-replicated.gro gmx-large.gro

#Center the bilayer in the box
editconf -f gmx-large.gro -o gmx-centered.gro -c yes

###___Solvate___###

#Output topology
python MDATopBuilder.py gmx-centered.gro  py-large-top.top

#Make a new topology with a position restraint file
python MDAPosTopBuilder.py gmx-centered.gro py-posre-top.top

#Energy minimize
grompp -f em-min.mdp -c gmx-centered.gro -p py-large-top.top -maxwarn 10 -o em-large.tpr
mdrun -deffnm em-large -v

#Position restraint minimization
grompp -f md-posre.mdp -c gmx-centered.gro -p py-posre-top.top -maxwarn 10 -o im-posre.tpr
mdrun -deffnm im-posre -v

#Production run
grompp -f md-prod.mdp -c im-posre.gro -p py-large-top.top -maxwarn 10 -o pr.tpr
mdrun -deffnm pr -v
