import sys

# infile=(open(sys.argv[1])).read()
infile = sys.argv[1]
outfile = sys.argv[2]

#Import Data
# infile = 'dppc_bilayer.pdb'
# outfile = "dppc_bilayer.csv"


file = open(infile)
data = file.read()

#Remove PDB headers, footers, and other information. Note that this makes our code
# need to be used in conjunction with Composition_Change.py so that box sizes, titles, etc
# are consistent throughout. Since we are, this is no issue. We don't need to store these
# headers, etc. for now. But note that if we want to combine code, it would be good
# to store all these values and print them directly back into the final output file.

b = ["TER","ENDMDL","TITLE     MIXED BILAYER","TITLE     DPPC BILAYER",
     "REMARK    THIS IS A SIMULATION BOX",
     "CRYST1   63.191   64.610  100.548  90.00  90.00  90.00 P 1           1",
     "MODEL        1"]

for i in range(0,len(b)):
    data = data.replace(b[i],"")



#Replace all spaces with commas.
data = data.replace(" ", ",")

#Iterate to remove all repeated commas
data = data.replace(",,",",")
data = data.replace(",,",",")
data = data.replace(",,",",")
data = data.replace(",,",",")
data = data.replace(",,",",")
data = data.replace(",,",",")

#Remove trailing commas because these make the last item "FALSE"
data = data.replace(",0.00,",",0.00")

#Remove blank lines
data = "".join([s for s in data.strip().splitlines(True) if s.strip()])

#Export to CSV (Actually just a text file)
out = open(outfile,'w') # Create an output file and print our annotations/data
out.write(data)
out.close()