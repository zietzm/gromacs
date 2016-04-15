It just prints the same thing for different indices. Just get an ndx file for all
headgroups and we'll be good Posre123.itp is all NC3.

We first make an index file that selects the atoms we want to be position
restrained.
```
make_ndx -f dppc_bilayer.gro -o index.ndx << a NC3
```
then a NC3
then q

Then we generate a position restraint file.

genrestr -f dppc_bilayer.gro -n index.ndx -o posre.itp

Finally, we just include this file in our topology file

We should figure out how to automatically make an index file using a bash
bash script.
