# Generating a Lipid Bicelle

We are trying to simulate a lipid bicelle using Gromacs. In order to do so, we
will need to use a large number of lipids, on the order of 1000. Doing an
all-atom simulation becomes computationally expensive so we will use a coarse-grained
forcefield. The Martini forcefield maps 4:1 atoms to coarse-grained beads,
which allows us to sacrifice a small amount of detail for a large improvement in
compuational time.

We first pick the necessary lipids to generate a bicelle. Per Jiang, et al.,
bicelles are composed of a lipid mixture of short- and long-tailed lipids.
We will use DPPC and DBPC (per Jiang) in this simulation.

A simple way to get a mixture of two lipids in a bilayer is to edit a DPPC
bilayer to have a certain composition of DBPC. We can start with a 128 DPPC
bilayer from the Martini website, chop off the ends of several lipids, replicate
the bilayer, and energy minimize the system.

## Generating Mixed bilayer
First, we download dppc_bilayer.gro from Martini. Then, we start Gromacs and
convert dppc_bilayer.gro to a .pdb file for easier editing.
```
source /usr/local/gromacs/bin/GMXRC
editconf -f dppc_bilayer.gro -o dppc_bilayer.pdb
```
A tool I used initially that made things easier was to then convert the PDB file
into a csv file with no spaces. This ensured that I could easily import the
bilayer into python. In order to convert to csv, I simply opened the pdb in a
spreadsheet and saved it as a csv.

We use a Python script to edit a desired number of DPPC lipids into DBPC.
```python
import numpy as np
import csv

# We need 35% DBP. Since we already have 128 DPP, we will convert 44.8 ~ 45 or 35% to DBP.
number_to_be_converted = 45
infile = 'dppc_bilayer.csv'

data = np.genfromtxt(infile, delimiter=",", dtype=None)
ldat = data.size
counter0 = 0
dconc = number_to_be_converted*6 #We remove six beads to convert each DPPC to DBPC

#Remove the desired number of DPPC-unique beads to make DBPC.
for i in range(0,ldat):
    for j in ['C2A','C3A','C4A','C2B','C3B','C4B']:
        if counter0 < dconc:
            if (j in data[i]) == True:
                counter0 = counter0+1
                data = np.delete(data, (i))

#Change the edited rows of each molecule to read, DBP instead of DPP.
for k in range(0,dconc):
    if ('NC3' in data[k]) == True:
        for l in range(k,k+6):
            ((data[l])[3]) = 'DBP'

ldat = data.size # Update the size to reflect edits made.

# Print Output in PDB format. We could edit title, box size, etc.
ann=("TITLE     MIXED BILAYER" + '\n'
     "REMARK    THIS IS A SIMULATION BOX" +'\n'
     "CRYST1   63.191   64.610  100.548  90.00  90.00  90.00 P 1           1" +'\n'
     "MODEL        1" +'\n')
fmt="%0s%7s%5s%4s%6s%12s%8s%8s%6s%6s" # Print the tuples to be exactly spaced as pdb.
x=''
for i in range(0,ldat):
   x = x + (fmt % (tuple(data[i]))) + '\n' #We should print automatically to PDB

out = open("Out.pdb",'w') # Create an output file and print our annotations/data
out.write(ann)
out.write(x)
out.close()

#Generate topology file
counter1 = 0
counter2 = 0
counter3 = 0
for m in range(0,ldat):
    if ('DPP' in data[m]) == True:
        counter1 = counter1 + 1
    if ('DBP' in data[m]) == True:
        counter2 = counter2 + 1
    if ('W' in data[m]) == True:
        counter3 = counter3 +1
ndp = str(counter1/12) # Account for the number of beads in each molecule.
ndb = str(counter2/6)
nwat = str(counter3)

top = open('topol.top', 'w')
top.write('#include "martini_v2.1.itp"'+ '\n'
           '#include "dppc_single.itp"'+'\n'
           '#include "dbpc_single.itp"'+'\n'
           '\n'
          '[ system ]' +'\n'
          'LIPID BICELLE SELF-ASSEMBLY'+'\n'
          '\n'
          '[ molecules ]'+'\n'
          'DPPC %s' % ndp +'\n'
          'DBPC %s' % ndb +'\n'
          'W %s' % nwat +'\n')
top.close()
```
The above code took our 128 DPPC bilayer, made a determined number of edits, and
generated a topology file with the exact composition of the output file.

To make numbering consistent throughout, we can use Gromacs to convert our
edited file back to .gro. `editconf -f Out.pdb -o Out.gro` Visually inspecting
the file, we notice that all numbering is consistent and the mixture is
as-specified.

We now want to increase the box size, increase the number of lipids, but keep
our initial DPPC/DBPC mixture ratio.

Just to be safe, minimize the box
```
;genbox -cp Out.gro -cs water.gro -o Out1.gro -maxsol 2000 -vdwd 0.21 ;Only if W was removed.
grompp -f minimization.mdp -c Out1.gro -p topol.top -maxwarn 10 -o Out1_min.tpr
mdrun -deffnm Out1_min -v
```
For our replication to proceed without error, I noticed that it is advantageous
to remove all water. We can do this manually with a text editor. Once it is done,
we need to manually edit the topology to reflect the lack of water. Rename the
waterless box to be Out2.gro.

Replicate the box nine times total. 3 in the x and 3 in the y directions. We
do not replicate in the z direction.
`genconf -f Out2.gro -o Out3.gro -nbox 3 3 1`

With a new concentration, we must edit the topology file. Simply multiply the
previous numbers by nine each.

At this point, we want to manually edit the .gro file to increase box size in the
x and y directions. Since we don't have any water at this point, we can also
reduce our box size in the z direction, so long as we don't interfere with the lipids.
Our new box will be called Out4.gro

For minimizations, we need water in the box. Simply solvate the box With 18,000
water-- the same concentration as the Martini file.
```
genbox -cp Out4.gro -cs water.gro -maxsol 18000 -vdwd 0.21 -o Out5.gro
grompp -f minimization.mdp -c Out5.gro -p topol.top -maxwarn 10 -o Out5_min.tpr
mdrun -deffnm Out5_min -v -nt 1
```
Now for our production run. This will take a considerable amount of time, on
the order of several hours.
```
grompp -f martini_md.mdp -c Out5_min.gro -p topol.top -maxwarn 10 -o Out5_Martini.tpr
mdrun -deffnm Out5_Martini -v
```
;;;; Compute densities across the box
;;; g_density -f Mixed5_Martini.gro -s Mixed5_Martini.tpr -o Mixed5_martini_Density.xvg
