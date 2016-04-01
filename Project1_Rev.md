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
We use a Python script to edit a desired number of DPPC lipids into DBPC.
```
import numpy as np
import csv

# We need 35% DBP. Since we already have 128 DPP, we will convert 44.8 ~ 45 or 35% to DBP.
number_to_be_converted = 45
infile = 'CM1.pdb'

data = np.genfromtxt(infile, delimiter=",", dtype=None)
ldat = data.size
counter0 = 0
dconc = number_to_be_converted*6 #We remove six beads to convert each DPPC to DBPC

for i in range(0,ldat):
    for j in ['C2A','C3A','C4A','C2B','C3B','C4B']:
        if counter0 < dconc:
            if (j in data[i]) == True:
                counter0 = counter0+1
                data = np.delete(data, (i))

#Edit the remaining rows of each molecule to read, DBP
for k in range(0,dconc):
    if ('NC3' in data[k]) == True:
        for l in range(k,k+6):
            ((data[l])[3]) = 'DBP'


ldat = data.size
# Print Output in PDB format
ann=("TITLE     MIXED BILAYER" + '\n'
     "REMARK    THIS IS A SIMULATION BOX" +'\n'
     "CRYST1   63.191   64.610  100.548  90.00  90.00  90.00 P 1           1" +'\n'
     "MODEL        1" +'\n')
fmt="%0s%7s%5s%4s%6s%12s%8s%8s%6s%6s"
x=''
for i in range(0,ldat):
   x = x + (fmt % (tuple(data[i]))) + '\n' #We should print automatically to PDB

out = open("OutT0.pdb",'w')
out.write(ann)
out.write(x)
out.close()

## csvfile = "</home/michael/Python/"
## with open('Cust_Martini.csv', "w") as output:
##     cwriter = csv.writer(output, delimiter=',', lineterminator='\n')
##     cwriter.writerows(data)

# We need to produce topology, as well.
# Update variable ldat since we removed many values


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
ndp = str(counter1/12)
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

;;;
;;;
;;;dppc_Martini.gro is the same as Martini's "dppc_bilayer.gro" Directly from Martini.

source /usr/local/gromacs/bin/GMXRC
; Convert bilayer from Martini into a pdb, then edit it to have desired composition.
; For ease, we will probably have removed waters in the pdb file.
; Convert the manually-created mixed bilayer from pdb to gro
editconf -f dppc_Martini.pdb -o Mixed0.gro

; Just to be safe, minimize the box
genbox -cp Mixed0.gro -cs water.gro -o Mixed1.gro -maxsol 2000 -vdwd 0.21
grompp -f minimization.mdp -c Mixed1.gro -p topol.top -maxwarn 10 -o Mixed1_min.tpr
mdrun -deffnm Mixed1_min -v

; Remove waters and edit topology (Mixed2.gro)

; Replicate the box DONE
genconf -f Mixed2.gro -o Mixed3.gro -nbox 3 3 1

; Edit topology DONE

; Increase box size MANUALLY -o Mixed4.gro

; Solvate
genbox -cp Mixed4.gro -cs water.gro -maxsol 18000 -vdwd 0.21 -o Mixed5.gro
grompp -f minimization.mdp -c Mixed5.gro -p topol.top -maxwarn 10 -o Mixed5_min.tpr
mdrun -deffnm Mixed5_min -v -nt 1

; Production run
grompp -f martini_md.mdp -c Mixed5_min.gro -p topol.top -maxwarn 10 -o Mixed5_Martini.tpr
mdrun -deffnm Mixed5_Martini -v

; Compute densities across the box
g_density -f Mixed5_Martini.gro -s Mixed5_Martini.tpr -o Mixed5_martini_Density.xvg

; Extend run a bit
grompp -f martini_md1.mdp -c Mixed5_Martini.gro -p topol.top -maxwarn 10 -o Mixed5_Martini2.tpr
mdrun -deffnm Mixed5_Martini2 -v

;;; Solvate CURRENT
;;; genbox -cp Mixed3.gro -cs water.gro -o Mixed4.gro -maxsol 18000 -vdwd 0.21
;;; Edit topology
;;; Minimization WILL NOT RUN
;;; grompp -f minimization.mdp -c Mixed4.gro -p topol.top -maxwarn 10 -o Mixed4_min.tpr
;;; mdrun -deffnm Mixed4_min -v
;;; Remove water and edit topology
