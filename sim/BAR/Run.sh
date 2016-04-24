
source /usr/local/gromacs/bin/GMXRC
#Solvate
genbox -cp gmx-nowater.gro -cs str-water.gro -vdwd 0.21 -o gmx-full-water.gro
