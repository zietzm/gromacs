source /usr/local/gromacs/bin/GMXRC

#Change composition and output topology
python py-comp-change.py str-dppc-bilayer.gro gmx-mixed-bilayer
python py-top-builder.py gmx-mixed-bilayer.gro gmx-mixed-top.top

#Energy minimize
grompp -f em.mdp -c gmx-mixed-bilayer.gro -p gmx-mixed-top.top -maxwarn 10 -o em-mix-bil.tpr
mdrun -deffnm em-mix-bil -v

#Replicate bilayer
genconf -f em-mix-bil.gro -o gmx-replicated.gro -nbox 3 3 1

#Change box size
python py-change-box.py gmx-replicated.gro gmx-large.gro

#Center the bilayer in the box
editconf -f gmx-large.gro -o gmx-centered.gro -center 12.5 12.5 5.25

#Remove water and then add it evenly.
python py-remove-water.py gmx-centered.gro gmx-nowater.gro
genbox -cp gmx-nowater.gro -cs str-water.gro -vdwd 0.21 -o gmx-full-water.gro

#Reorder gmx-full-water.gro so that it works for the topology
python py-reorder.py gmx-full-water.gro gmx-ordered

#Output topology
python py-top-builder.py gmx-full-water.gro gmx-large-top.top

#Make position restrained topologies
python py-posre-top-builder.py gmx-full-water.gro gmx-posre-top.top
python py-posre-all-top-builder.py gmx-full-water.gro gmx-posre-all-top.top

#Restrain all lipid atoms and let the water molecules readjust.
grompp -f em-max.mdp -c gmx-ordered.gro -p gmx-posre-all-top.top -maxwarn 10 -o em-watermin.tpr
mdrun -deffnm em-watermin -v #-nt 1

#Energy minimize
grompp -f em-max.mdp -c em-watermin.gro -p gmx-posre-top.top -maxwarn 10 -o em-prm.tpr
mdrun -deffnm em-prm -v #nt 4 #Blowing up. Need nt 1

#Position restraint minimization
grompp -f md-posre.mdp -c em-prm.gro -p gmx-posre-top.top -maxwarn 10 -o em-posre.tpr
mdrun -deffnm em-posre -v

#Production run
grompp -f md-martini.mdp -c em-posre.gro -p gmx-large-top.top -maxwarn 10 -o md-pr.tpr
mdrun -deffnm md-pr -v