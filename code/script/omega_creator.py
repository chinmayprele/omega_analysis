"""
- finds the omega dN/dS for a sliding window

- INPUT:
    - FASTA of dNTPS (from directory)
    - request window size
- OUTPUT:
    - file of oemga values for each sliding window


FINISH CREATING THIS file

READ THE CODEML FROM:
    ~/Downloads/paml4.8/bin/codeml
IT READS IN THE codeml.ctl AT:
    ~/Downloads/paml4.8/codeml.ctl
NEED TO OUTPUT INTO A FORMAT THAT IT CAN READ.
PARAMETERS ALREADY SET; DO NOT CHANGE

"""

### imports ==================================================================
from Bio import SeqIO                   # to read fasta in and out
from pprint import pprint as pprint     # testing
import sys                              # script creation
import os                               # to run codeml from command line
import subprocess
import time
from tqdm import tqdm

### inputs ==================================================================
# input_file = input( "input ALN file path: " )
# input_file = "../data/test/ilp1_fna.aln"
# window_size = int(input( "window size: " ))
window_size = 10            # in terms of codons

codeml_path = "~/Downloads/paml4.8/bin/codeml"

### code ==================================================================

codon_dict = {}
aa_length = 0

with open( "../../data/tor/Tor.codons", "r") as infile:
    for line in tqdm(infile, desc="Reading in codons file", ascii=True):
        line = line.split()
        if line[0] not in codon_dict.keys():
            codon_dict[ line[0] ] = line[1:]
            aa_length = len( codon_dict[ line[0] ] )

# FILTER "???"
# for item in codon_dict.keys():
#     codon_dict[item].remove( "???" )

omega_lst = [[
    "range_start",
    "range_stop",
    "omega",
    "kappa",
    "window_size"
]]

highly_evolvable = {}

with open( "omegas.lst", "w" ) as omfile:
    print( "\t".join( omega_lst[0] ), file=omfile )
    time.sleep(1)
    for i in tqdm(range( 0, aa_length-window_size+1 ), desc="Calculate omega", ascii=True):
    # for i in range( 0, aa_length-window_size+1 ):
        start = i
        stop = i+window_size
        temp_comparisons = []
        with open( "example.aln", "w" ) as outfile:
            for item in codon_dict.keys():
                print( ">{0}\n{1}".format( item, "".join(codon_dict[item][ start:stop ]) ), file=outfile )
                temp_comparisons.append( ">{0}\n{1}".format( item, "".join(codon_dict[item][ start:stop ]) ) )
        os.system( "{0} {1}".format( codeml_path, ">/dev/null 2>&1" ) )
        try:
            direct_output = subprocess.check_output("cat Tor_mlc | grep \"omega\"", shell=True)
            omega = direct_output.strip().split()[-1].decode("utf-8")
            # print( "="*10,start,"="*10, stop,"="*10, omega  )
            # omega_lst.append( [start+1, stop+1, omega] )
            direct_output = subprocess.check_output("cat Tor_mlc | grep \"kappa\"", shell=True)
            kappa = direct_output.strip().split()[-1].decode("utf-8")
            # print( "="*10,start,"="*10, stop,"="*10, kappa  )
            print( "\t".join([str(start+1), str(stop+1), str(omega), str(kappa), str(window_size)]), file=omfile )
            if float(omega) > 1:
                highly_evolvable[ "{0}:{1}".format( start, stop ) ] = temp_comparisons
            time.sleep(1)
            # print("here")
        except subprocess.CalledProcessError as e:
            # raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
            continue
        # sys.exit( "creation" )

with open( "Tor_high_omegas.fna", "w" ) as outfile:
    for item in highly_evolvable.keys():
        print( ">REGION{0}\n{1}".format( item, highly_evolvable[item], file=outfile, end="\n\n" ) )
