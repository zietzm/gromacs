source /usr/local/gromacs/bin/GMXRC

python martinize.py -f 4i1q.pdb -x str-martini-NBAR.pdb -o atop.top
python insane.py -f str-martini-NBAR.pdb -o gmx-BAR-lipids.gro -p topol.top -l DPPC:3 -l DTPC:1 -u DPPC:3 -u DTPC:1 -pbc square -box 25.0,25.0,10.0 -center -sol W

grompp -f em.mdp -c gmx-BAR-lipids.gro -p topol.top -maxwarn 10 -o em-prolip.tpr
mdrun -deffnm em-prolip -v

grompp -f md-martini.mdp -c em-prolip.gro -p topol.top -maxwarn 10 -o md-prolip.tpr
mdrun -deffnm md-prolip -v


