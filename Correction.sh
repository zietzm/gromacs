#!/bin/bash
source /usr/local/gromacs/bin/GMXRC

genbox -cp vmBicelle11Big.gro -cs water.gro -maxsol 15000 -vdwd 0.21 -o v11Solv.gro

grompp -f Posremin.mdp -c v11Solv.gro -p Solv.top -maxwarn 10 -o v11pr.tpr
mdrun -deffnm v11pr -v

grompp -f martini_extend.mdp -c v11pr.gro -p SolvNoPr.top -maxwarn 10 -o v11Martini.tpr
tmux new-session -d -s Martini-Run 'mdrun -deffnm v11Martini -v'
