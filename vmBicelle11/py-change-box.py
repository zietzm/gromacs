import sys

infile = sys.argv[1]
outfile = sys.argv[2]

# infile = "Out1.gro"
# outfile = "Out2.gro"

file=open(infile)
data=file.read()

data = data.replace("18.95730  19.38300  10.05480","25.00000   25.00000  10.00000") #Make this a general removal instead of specific

out = open(outfile,'w') # Create an output file and print data
out.write(data)
out.close()
