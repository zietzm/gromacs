
#!/bin/bash
ip=dppc_bilayer.csv
op=OutTest.pdb
oc=OutTest.csv

python Composition_Change.py $ip $op $oc
