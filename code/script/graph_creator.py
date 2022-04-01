"""
- creates graph from omegas.lst
- Outputs another file with sequence of regions with high omega values (>1)
"""

### imports ==================================================================
from pprint import pprint as pprint
import matplotlib.pyplot as plt
import numpy as np
import statistics
import subprocess
from Bio import SeqIO
import sys

### variables ==================================================================

gene = "Ilp7"
omegas_file = "/Users/rele.c/Desktop/github_repos/omega_analysis/model_data/{0}/omegas.lst".format( gene )

out_image_path = "/Users/rele.c/Desktop/github_repos/omega_analysis/model_data/{0}/images".format( gene )

# codon_file = "../data/chico/chico.codons"
# peptide_file = "../data/chico/chico.pep"

omega_cutoff = 1

# high_omega_file = "_{0}_highomegas_{1}.seq".format( gene, str(omega_cutoff) )


### code ==================================================================

# start, stop, omega, kappa, window_size
omega_lst = []

with open( omegas_file,  "r" ) as infile:
    for line in infile:
        if "range_start" in line:
            continue
        omega_lst.append( [float(i) for i in line.strip().split()]  )

### create graph ==================================================================
starts, omegas = [], []
for item in omega_lst:
    starts.append( item[0] )
    omegas.append( item[2] )

fig = plt.figure(figsize=(30,3))
axes = fig.add_axes([0.1,0.1,0.8,0.8])
x = starts
y = omegas
axes.plot(x,y)

axes.set_ylim([0,5])

plt.title( gene )
plt.grid(color='k', linestyle='-', linewidth=0.1)
plt.xlabel( "peptide position" )
plt.ylabel( "dN/dS" )

# plt.show()
plt.savefig( "{0}/dNdS_graph_{1}.png".format( out_image_path, gene ) )

sys.exit()

### extract sequences ==================================================================

with open( "temp.tab", "w" ) as tmpfile:
    for line in omega_lst:
        if line[2] > omega_cutoff:
            tmp = [ "chr" ]
            tmp.append( str(int(line[0])) )
            tmp.append( str(int(line[1])) )
            tmp.append( str(float(line[2])) )
            tmp.append( str(float(line[3])) )
            tmp.append( str(int(line[4])) )
            print( "\t".join(tmp), file=tmpfile )

out = subprocess.Popen(['mergeBed', '-i', 'temp.tab', '-c', '4,5,6', '-o', 'mean,mean,distinct'],
           stdout=subprocess.PIPE,
           stderr=subprocess.STDOUT)
stdout,stderr = out.communicate()
lines = stdout.decode("utf-8").strip().split("\n")

ranges = []

for item in lines:
    t = []
    t.append( int(item.split()[1]) )
    t.append( int(item.split()[2]) )
    ranges.append( t )

with open( "codon"+high_omega_file, "w" ) as outfile:
    with open( codon_file, "r" ) as codfile:
        for item1 in codfile:
            temp = item1.split()[1:]
            for item in ranges:
                print( "\t".join([str(item[0]), str(item[1]), "".join(temp[ item[0]: item[1] ])]), file=outfile )
            break

with open( "pep"+high_omega_file, "w" ) as outfile:
    records = list(SeqIO.parse( peptide_file , "fasta"))
    for item in ranges:
        print( "\t".join([str(item[0]), str(item[1]), "".join(records[0].seq[ item[0]: item[1] ])]), file=outfile )
