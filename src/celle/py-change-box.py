import sys

infile = sys.argv[1]
outfile = sys.argv[2]

file=open(infile)
data=file.read()

data = data.replace(data[-29:],"25.00000  25.00000  12.50000")

out = open(outfile,'w') 
out.write(data)
out.close()
