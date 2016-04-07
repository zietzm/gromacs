import sys

infile = sys.argv[1]
outfile = sys.argv[2]

# infile = "Out1.gro"
# outfile = "Out2.gro"

file=open(infile)
data=file.read()

data = data.replace("6.31910   6.46100  10.05480","23.90000   20.40000  8.10000")

out = open(outfile,'w') # Create an output file and print data
out.write(data)
out.close()
