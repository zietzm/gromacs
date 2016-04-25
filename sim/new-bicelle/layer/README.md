Required files to generate a bilayer of DXPC:

* bash-gen-bilayer.sh
* em.mdp
* md-bilayer.mdp
* py-top-builder.py
* str-dxpc-single.gro
* str-water.gro
* top-dtpc-single.itp
* top-dxpc-single.itp
* top-martini-v2.1.itp

Output files:
* gmx-bilayer.gro
* gmx-bilayer.top

Generating this bilayer should not take more than ten minutes.

I recommend visually inspecting the bilayer with VMD before moving on. Sometimes, the short simulation
time will lead to a deformed bilayer, etc. By keeping the simulation short, this becomes a possibility.
Simply restart the simulation from the beginning, and you should get a workable bilayer.

Below, annotated bash file:
```
# Add 128 DXPC
genbox -ci str-dxpc-single.gro -box 6.31915 6.46099 10.05482 -nmol 128 -try 1000 -o gmx-vac-lipids.gro

#Generate topology
python py-top-builder.py gmx-vac-lipids.gro gmx-vac.top

#Energy minimize
grompp -f em.mdp -c gmx-vac-lipids.gro -p gmx-vac.top -maxwarn 10 -o em-vac.tpr
mdrun -deffnm em-vac -v

#Solvate and increase box size to be enough for 2000 lipids
genbox -cp em-vac.gro -cs str-water.gro -vdwd 0.21 -maxsol 2000 -box 6.5 6.5 12.5 -o gmx-lipids.gro

#Generate topology
python py-top-builder.py gmx-lipids.gro gmx-lipids.top

#Minimization
grompp -f em.mdp -c gmx-lipids.gro -p gmx-lipids.top -maxwarn 10 -o em-lipids.tpr
mdrun -deffnm em-lipids -v

#Production run
grompp -f md-bilayer.mdp -c em-lipids.gro -p gmx-lipids.top -maxwarn 10 -o md-lipids.tpr
mdrun -deffnm md-lipids -v

#Put all molecules into the frame and output fixed coordinate file.
echo "0" | trjconv -f md-lipids.xtc -s md-lipids.tpr -pbc mol -o gmx-bilayer.xtc
echo "0" | trjconv -f gmx-bilayer.xtc -s md-lipids.tpr -dump 30000 -vel -o gmx-bilayer.gro

cp gmx-lipids.top gmx-bilayer.top

# #Delete Unnecessary files and reorganize
# rm em-lipids.edr em-lipids.gro em-lipids.log em-lipids.tpr tm-lipids.trr gmx-lipids.gro gmx-lipids.xtc
# rm md-lipids.cpt md-lipids.edr md-lipids.gro md-lipids.log md-lipids.tpr md-lipids.xtc mdout.mdp
```
