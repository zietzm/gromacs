#!/bin/bash
set -e

source /usr/local/gromacs/bin/GMXRC

python py-change-comp.py str-dppc-bilayer.gro gmx-mixed-bilayer.gro
python py-gen-top.py gmx-mixed-bilayer.gro gmx-mixed-bilayer.top

grompp -f em.mdp ####
mdrun -deffnm em-mix-bil -v

genconf -f em-mix-bil.gro -o gmx-replicated.gro -nbox 3 3 1
python py-change-box.py gmx-replicated.gro gmx-large.gro
editconf -f gmx-large.gro -o gmx-centered.gro -center 12.5 12.5 6.25
python py-remove-water.py gmx-centered.gro gmx-nowater.gro
genbox -cp gmx-nowater.gro -cs str-water.gro -vdwd 0.21 -o gmx-full-water.gro

grompp -f em.mdp -c gmx-ordered.gro -p gmx-posre-all-top.top -maxwarn 10 -o em-watermin.tpr
mdrun -deffnm em-watermin -v 
grompp -f em.mdp -c em-watermin.gro -p gmx-posre-top.top -maxwarn 10 -o em-prm.tpr
mdrun -deffnm em-prm -v
grompp -f md-posre.mdp -c em-prm.gro -p gmx-posre-top.top -maxwarn 10 -o em-posre.tpr
mdrun -deffnm em-posre -v
grompp -f md-martini.mdp -c em-posre.gro -p gmx-large-top.top -maxwarn 10 -o md-pr.tpr
tmux new-session -d -s prod_run 'mdrun -deffnm md-pr -v'
