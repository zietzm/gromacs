#!/bin/bash
source /usr/local/gromacs/bin/GMXRC
make_ndx -f dppc_bilayer.gro -o index.ndx
a NC3
q

