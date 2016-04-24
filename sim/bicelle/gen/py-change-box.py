import sys

infile = sys.argv[1]
outfile = sys.argv[2]

file=open(infile)
data=file.read()

data = data.replace(data[-29:],"25.00000  25.00000  10.50000")

out = open(outfile,'w') # Create an output file and print data
out.write(data)
out.close()
