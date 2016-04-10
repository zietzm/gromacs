# Using MDAnalysis
MDAnalysis (mda) will make our simulation so much simpler than I had previously thought possible. Mda will facilitate mainpulation 
of files without constantly changing between gro, pdb, and csv. I have come up with the following outline of a bicelle generation
using mda instead of several disconnected python scripts. I should note, however, that while my python scripts were far from 
perfect, they were-- for the most part-- working quite well for the tasks prescribed.

1. Change composition of dppc-bilayer.gro either while it is gro or after converting to pdb.
2. Run energy minimization
3. Remove water
4. Replicate bilayer
5. Change box size
6. Solvate
7. Run energy minimization
8. Production run

This new tool will-- at absolute minimum at this moment-- facilitate atom selection in various file formats. At very minimum, this
capability makes the most tricky python script obsolete, PDB_to_CSV.py. 