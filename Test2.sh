#!/bin/bash

# ip=dppc_bilayer.csv
# op=dppc_bilayer.pdb
# oc=OutTest.csv
# og=test123.gro
#
# source /usr/local/gromacs/bin/GMXRC

tmux new-session -d -s my_session 'cmatrix'
tmux detach -s my_session

#tmux kill-session -t my_session
