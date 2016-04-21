#!/bin/bash

#Change composition
python py-comp-change.py str-dppc-bilayer.gro gmx-mixed-bilayer

#Output topology
python py-top-builder.py gmx-mixed-bilayer.gro gmx-mixed-top.top

#Energy minimize
grompp -f em.mdp -c gmx-mixed-bilayer.gro -p gmx-mixed-top.top -maxwarn 10 -o em-mix-bil.tpr
mdrun -deffnm em-mix-bil -v 

#Replicate bilayer
genconf -f em-mix-bil.gro -o gmx-replicated.gro -nbox 3 3 1

#Change box size
python py-change-box.py gmx-replicated.gro gmx-large.gro

#Center the bilayer in the box
editconf -f gmx-large.gro -o gmx-centered.gro -c yes

#Output topology
python py-top-builder.py gmx-centered.gro gmx-large-top.top

#Make a new topology with a position restraint file
python py-posre-top-builder.py gmx-centered.gro gmx-posre-top.top

#Energy minimize
grompp -f em.mdp -c gmx-centered.gro -p gmx-large-top.top -maxwarn 10 -o em-large.tpr
mdrun -deffnm em-large -v 
    
#Position restraint minimization
grompp -f md-posre.mdp -c gmx-centered.gro -p gmx-posre-top.top -maxwarn 10 -o em-posre.tpr
mdrun -deffnm em-posre -v
    
#Production run
grompp -f md-prod.mdp -c em-posre.gro -p gmx-large-top.top -maxwarn 10 -o md-pr.tpr
mdrun -deffnm md-pr -v

