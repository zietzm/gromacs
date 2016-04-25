#!/bin/bash
set -e

### Requires the following files:
### bash-gen-bicelle.sh
### em.mdp
### md-martini.mdp
### md-posre.mdp
### py-change-box.py
### py-comp-change.py
### py-posre-all-top-builder.py
### py-posre-top-builder.py
### py-remove-water.py
### py-reorder.py
### py-top-builder.py
### README.md
### str-dxpc-bilayer.gro (copy of gmx-bilayer.gro)
### str-water.gro
### top-dtpc-posre.itp
### top-dtpc-posre-all.itp
### top-dtpc-single.itp
### top-dxpc-posre.itp
### top-dxpc-posre-all.itp
### top-dxpc-single.itp
### top-martini-v2.1.itp

source /usr/local/gromacs/bin/GMXRC
python py-comp-change.py str-dxpc-bilayer.gro gmx-mixed-bilayer
python py-top-builder.py gmx-mixed-bilayer.gro gmx-mixed-top.top
grompp -f em.mdp -c gmx-mixed-bilayer.gro -p gmx-mixed-top.top -maxwarn 10 -o em-mix-bil.tpr
mdrun -deffnm em-mix-bil -v -nt 1
genconf -f em-mix-bil.gro -o gmx-replicated.gro -nbox 3 3 1
python py-change-box.py gmx-replicated.gro gmx-large.gro
editconf -f gmx-large.gro -o gmx-centered.gro -center 12.5 12.5 6.25
python py-remove-water.py gmx-centered.gro gmx-nowater.gro
genbox -cp gmx-nowater.gro -cs str-water.gro -vdwd 0.21 -o gmx-full-water.gro
python py-reorder.py gmx-full-water.gro gmx-ordered
python py-top-builder.py gmx-full-water.gro gmx-large-top.top
python py-posre-top-builder.py gmx-full-water.gro gmx-posre-top.top
python py-posre-all-top-builder.py gmx-full-water.gro gmx-posre-all-top.top
grompp -f em.mdp -c gmx-ordered.gro -p gmx-posre-all-top.top -maxwarn 10 -o em-watermin.tpr
mdrun -deffnm em-watermin -v 
grompp -f em.mdp -c em-watermin.gro -p gmx-posre-top.top -maxwarn 10 -o em-prm.tpr
mdrun -deffnm em-prm -v
grompp -f md-posre.mdp -c em-prm.gro -p gmx-posre-top.top -maxwarn 10 -o em-posre.tpr
mdrun -deffnm em-posre -v
grompp -f md-martini.mdp -c em-posre.gro -p gmx-large-top.top -maxwarn 10 -o md-pr.tpr
tmux new-session -d -s prod_run 'mdrun -deffnm md-pr -v'
