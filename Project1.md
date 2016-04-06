# Lipid Bicelle Simulation

We are trying to generate a bicelle using Gromacs. Because a large number of
lipids are required, computational time becomes arbitrarily large. To fix this
problem, we will use the Martini Coarse-Grained force field to simplify the
computation.

First, we need to pick two lipids to use. Bicelles need a both a short- and a
long-tailed lipid. Next, we need to pick the composition of our mixture: we want
a favorable mix of short- and long-tailed lipids to give us the best chance of
forming a bicelle. Once we have a composition of two lipids, we need to pick a
box size for our lipid simulation.

*Note: It is very possible that we received strange results in our first
simulation because we had a box that was too small with periodic boundary
conditions. We should try to make the box a better size this time.*

<img src ="/images/lipid-md2.png" alt = "Lipid Image" style="width:300px;height:150px;">

In 2007, Yong Jiang, in his Ph.D. dissertation at Emory University, showed how a
lipid bicelle could be created with various types of short- and long-tailed
lipids. Specifically, he used DPPC (16 carbon tails) and DBPC (4 carbon tails).
A bicelle was generate for a concentration of 65% DPPC (628 DPPC and 338 DBPC
in a box of size 23.9 8.1 20.4). This simulation used an integration time of 20
fs over a hundreds-of-ns simulation time. This gives us around 12,500,000
steps needed. Jiang used three-dimensional periodic boundary conditions. He used
a "shift function to coulombic force when the interaction distance falls between
0 nm and 1.2 nm" to make the electric field "smoothly" go to zero. Van der Waals
forces were treated the same way, but the shift function was, "on between 0.9 nm
and 1.2 nm." Jiang noted that DPPC behaves like a fluid above 315K, and he
chose the simulation to take place at 323K.

As in Jiang, we will use the steepest descent algorithm in "minimization.mdp" to
minimize the system before it is sent to the "production run".


If all the above have been chosen-- lipids, compositions, box size,
number-- we need to start using Gromacs.

## Running simulation remotely

Because any complex simulation will require significant computing power and time,
we do not generally run simulations on personal computers. In fact, simulations
will be run almost exclusively through remote supercomputers or high-powered
workstations.

Working remotely has many advantages and is not difficult. Gromacs is
especially easy to work remotely because the entire program runs in the linux
terminal. It is for this reason that running simulations remotely is easiest if
you are working with a linux machine to begin with.

If you are using a linux personal computer (I use Ubuntu):

1. Open a terminal
2. Start ssh with the following command, and, when prompted, enter your password.
```
~$ ssh USER@ADDRESS
```
3. Transfer files to the remote station. If both the local and remote system are
running linux, we can add our files using Nautilus, the system file manager.
Simply open Nautilus and login using "Connect to Server". This makes the remote
computer's file system as accessible as if it were a local storage device. Now
we simply copy and paste our files to the remote directory.
4. Go back to terminal and CD into the directory you put the files in step 3.
5. Make sure tmux is installed so that we don't have to keep our terminal and
SSH session open during the entire course of the simulation.

If you are using a Windows personal computer:

We will need several tools in order to have as seamless an experience as we would
on a linux machine. First, we need to download PuTTY-- an SSH client-- and WinSCP,
a remote file manager that will act as a substitute for Nautilus's built-in
functionality.

## Setup lipids-water system

First, open Gromacs. Under my installation we enter:
```
source /usr/local/gromacs/bin/GMXRC
```

Let's start by making a box and filling it with one of our chosen lipids.
```
genbox -ci dppc_single.gro -nmol 628 -box 23.9 8.1 20.4 -try 500 -o lipid1.gro
```
Above, -ci is the molecule as a .pdb or a .gro, -nmol is the number of
molecules, -box gives
the box dimensions in nanometers, -try is the number of attempts to
generate the box (sometimes trying to throw a large number of lipids
together into a box doesn't work), -o gives the output as a .gro file.


Now lets setup and run an equilibration to minimise the system's free
energy. <br><br> Here, -f is the md file, which tells mdrun what to do, -c is
input, -p is the topology file to be made, -maxwarn is optional to tell
Gromacs that we only want a certain number of warnings, and -o is, as
usual, our output.
```
grompp -f minimization.mdp -c lipid1.gro -p lipid.top -maxwarn 10 -o lipid1-min.tpr
mdrun -deffnm lipid1-min -v
```

That shouldn't have taken long at all. We now move on by adding our other lipid.
To avoid any complications that arise when we try to merge to output files, we
will use a trick to "solvate" the output file from our first simulation with our
chosen second lipid.
```
genbox -cp lipid1-min.gro -cs dbpc_single.gro -o lipids.gro -maxsol 338
```

Again, we setup and run a minimisation of the system's free energy. First,
though, we need to open our topology file, lipid.top and uncomment the lipid
we just added. You can do this with nano lipid.top or just open it using a text
editor.
```
grompp -f minimization.mdp -c lipids.gro -p lipid.top -maxwarn 10 -o both-min.tpr
mdrun -deffnm both-min -v
```

Now we have an equilibrated system with two types of lipids. To simulate real
conditions, we will add water to the system. This time, we solvate our last
output with water.
```
genbox -cp both-min.gro -cs water.gro -o water-lipids.gro -maxsol 768 -vdwd 0.21
```

Again, edit the topology file to uncomment water. Now run the equilibration.
```
grompp -f minimization.mdp -c water-lipids.gro -p lipid.top -maxwarn 10 -o water-lipids-min.tpr
mdrun -deffnm water-lipids-min -v
```

## Equilibrate lipids-water system

It is now time for us to run the full simulation using our real .mdp file, martini_md.mdp. This .mdp file will use 12,500,000 steps of 20 fs to simulate a time
of 250 ns.
```
grompp -f martini_md.mdp -c water-lipids-min.gro -p lipid.top -maxwarn 10 -o solv-md.tpr
tmux
mdrun -deffnm solv-md -v
```
Now we wait for the simulation to run. This could take a long time, but
luckily we ran the process in tmux. This allows us to press <kbd>Ctrl</kbd>+<kbd>b</kbd>
<kbd>d</kbd> to detach from the tmux session. To test that it worked, we can
type `tmux attach` to check that the process is running even when we close the
SSH session. If everything worked with tmux, we can again detach from tmux and
type `exit` to exit the SSH.

## Inspect product

We have now finished what should be our bicelle! To check how this simulation
went, lets use VMD to visually inspect our product.
```
vmd solv-md.gro
```

Lets draw the bonds so we can better see what's going on. Open the TK console
from the extensions tab and draw the bonds.
```
source cg_bonds.tcl
cg_bonds -top lipid.top
```
