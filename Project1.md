# Lipid Bicelle Simulationss

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

If all the above have been chosen-- lipids, compositions, box size,
number-- we need to start using Gromacs.

## Setup lipids-water system

First, open Gromacs. Under my installation we enter:
```
source /usr/local/gromacs/bin/GMXRC
```

Let's start by making a box and filling it with one of our chosen lipids.
<pre><code>genbox</code></pre>
Above, -ci is the molecule as a .pdb or a .gro, -nmol is the number of
molecules, -box gives
the box dimensions in nanometers, -try is the number of attempts to
generate the box (sometimes trying to throw a large number of lipids
together into a box doesn't work), -o gives the output as a .gro file.

<br><br>

Now lets setup and run an equilibration to minimise the system's free
energy. <br><br> Here, -f is the md file, which tells mdrun what to do, -c is
input, -p is the topology file to be made, -maxwarn is optional to tell
Gromacs that we only want a certain number of warnings, and -o is, as
usual, our output.
<pre><code>grompp
mdrun
</code></pre>

That shouldn't have taken long at all. We now move on by adding our other lipid.
To avoid any complications that arise when we try to merge to output files, we
will use a trick to "solvate" the output file from our first simulation with our
chosen second lipid.
<pre><code>genbox</code></pre>

Again, we setup and run a minimisation of the system's free energy. First,
though, we need to open our topology file, lipid.top and uncomment the lipid
we just added. You can do this with nano lipid.top or just open it using a text
editor.
<pre><code>grompp
mdrun</code></pre>

Now we have an equilibrated system with two types of lipids. To simulate real
conditions, we will add water to the system. This time, we use <pre><code> genbox </code></pre>
to solvate our last output with our system's actual solvent. Again, we run a
minimisation.
<pre><code>genbox</code></pre>

Again, edit the topology file to uncomment water. Now run the equilibration.
<pre><code>grompp
mdrun</code></pre>

## Equilibrate lipids-water system

It is now time for us to run the full simulation using our real .mdp file,
martini_md.mdp. This .mdp file will use XXXX steps of 30 fs to simulate a time
of 200 ns.
<pre><code>grompp
mdrun</code></pre>

## Inspect product

We have now finished what should be our bicelle! To check how this simulation
went, lets use VMD to visually inspect our product.
<pre><code>vmd </code></pre>

Lets draw the bonds so we can better see what's going on. Open the TK console
from the extensions tab and draw the bonds.
<pre><code>source cg_bonds.tcl
cg_bonds -top lipid.top</code></pre>
