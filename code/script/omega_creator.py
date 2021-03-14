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
import modules
from subprocess import Popen, PIPE

from Bio.Phylo.PAML import codeml


### inputs ==================================================================
# input_file = input( "input ALN file path: " )
# input_file = "../model_data/test/ilp1_fna.aln"
# window_size = int(input( "window size: " ))
window_size = 10            # in terms of codons
# window_size = input( "Window size: " )            # in terms of codons

codeml_path = "/Users/rele.c/Downloads/omega_analysis/paml4.8/bin/codeml"

### finding gene ==================================================================

from pathlib import Path

path = "/Users/rele.c/Downloads/omega_analysis/model_data/"

genes = []
for temp in Path(path).iterdir():
	genes.append( str(temp).split("/")[-1] )
genes.sort()

for i in range( len( genes ) ):
	print( "{0:>4}: {1}".format( i, genes[i] ) )

gene = genes[ int(input( "Which gene should be analyzed: " )) ]
# gene = genes[ 4 ]

aln_path = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/raw_data/{0}.msa.aln".format( gene )
codons_path = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/{0}.codons".format( gene )
mlc_path = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/{0}.mlc".format( gene )
tree_path = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/raw_data/{0}.nwk".format( gene )
gene_dir_path = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/".format( gene )
full_gene_analysis = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/full_gene_{0}.out".format( gene )
codeml_ctl_path = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/ctl_files/codeml.ctl".format( gene )
mini_ctl_path = "/Users/rele.c/Downloads/omega_analysis/model_data/{0}/ctl_files/mini_codeml.ctl".format( gene )

# tree_path = "/Users/rele.c/Downloads/omega_analysis/model_data/test/raw_data/test_transl_subalign.nwk"


print( "aln_path: ".rjust(20), aln_path )
print( "codons_path: ".rjust(20), codons_path )
print( "mlc_path: ".rjust(20), mlc_path )
print( "tree_path: ".rjust(20), tree_path )
print( "gene_dir_path: ".rjust(20), gene_dir_path )
print( "codeml_ctl_path: ".rjust(20), codeml_ctl_path )

modules.codon_creator( aln_path, codons_path )
time.sleep(1)

# sys.exit("test in_python codeml")

### code ==================================================================

codon_dict = {}
aa_length = 0

codeml_params = {
	"seqfile": aln_path,
	"treefile": tree_path,
	"outfile": mlc_path,
	"runmode": 0,
	"seqtype": 1,
	"CodonFreq": 2,
	"model": 0,
	"NSsites": 0,
	"fix_kappa": 0,
	"kappa": 1,
	"fix_omega": 0,
	"omega": 1,
	"ncatG": 10,
}

# create codeml.ctl file for entire gene
with open( codeml_ctl_path, "w" ) as ctl:
	print( "{0} = {1}".format( "seqfile", codeml_params["seqfile"] ), file=ctl )
	print( "{0} = {1}".format( "treefile", codeml_params["treefile"] ), file=ctl )
	print( "{0} = {1}".format( "outfile", codeml_params["outfile"] ), file=ctl )
	print( file=ctl )
	print( "{0} = {1}".format( "runmode", codeml_params["runmode"] ), file=ctl )
	print( "{0} = {1}".format( "seqtype", codeml_params["seqtype"] ), file=ctl )
	print( "{0} = {1}".format( "CodonFreq", codeml_params["CodonFreq"] ), file=ctl )
	print( "{0} = {1}".format( "model", codeml_params["model"] ), file=ctl )
	print( "{0} = {1}".format( "NSsites", codeml_params["NSsites"] ), file=ctl )
	print( "{0} = {1}".format( "fix_kappa", codeml_params["fix_kappa"] ), file=ctl )
	print( "{0} = {1}".format( "kappa", codeml_params["kappa"] ), file=ctl )
	print( "{0} = {1}".format( "fix_omega", codeml_params["fix_omega"] ), file=ctl )
	print( "{0} = {1}".format( "omega", codeml_params["omega"] ), file=ctl )
	print( "{0} = {1}".format( "ncatG", codeml_params["ncatG"] ), file=ctl )
time.sleep( 1 )

# codeml_ctl_path = "/Users/rele.c/Downloads/omega_analysis/model_data/test/bare_codeml.ctl"

# run_line = "{0} {1} >/dev/null 2>&1".format( codeml_path, codeml_ctl_path )
# print( "run_line: ", run_line )
# codeml_run = Popen([run_line], stdin=PIPE, shell=True)
# time.sleep(1)
# codeml_run.communicate( input = "\n".encode('utf-8') )

# add the omega and cappa to a file
with open( full_gene_analysis, "w" ) as outfile:
	direct_output = subprocess.check_output("cat {0} | grep \"omega (dN/dS)\"".format( mlc_path ), shell=True)
	omega = direct_output.strip().split()[-1].decode("utf-8")
	direct_output = subprocess.check_output("cat {0} | grep \"kappa (ts/tv)\"".format( mlc_path ), shell=True)
	kappa = direct_output.strip().split()[-1].decode("utf-8")
	print( "{0} omega (dN/dS) = {1}".format( gene, omega ), file=outfile )
	print( "{0} kappa (ts/tv) = {1}".format( gene, kappa ), file=outfile )
	print( "Added full gene data to {0}".format( full_gene_analysis ) )

with open( codons_path, "r") as infile:
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

# create temp ctl file for codeml

temp_mlc = gene_dir_path + "temp.mlc"
temp_aln = gene_dir_path + "temp.aln"

codeml_params = {
	"seqfile": temp_aln,
	"treefile": tree_path,
	"outfile": temp_mlc,
	"runmode": 0,
	"seqtype": 1,
	"CodonFreq": 2,
	"model": 0,
	"NSsites": 0,
	"fix_kappa": 0,
	"kappa": 1,
	"fix_omega": 0,
	"omega": 1,
	"ncatG": 10,
}

# create codeml.ctl file for entire gene
with open( mini_ctl_path, "w" ) as ctl:
	print( "{0} = {1}".format( "seqfile", codeml_params["seqfile"] ), file=ctl )
	print( "{0} = {1}".format( "treefile", codeml_params["treefile"] ), file=ctl )
	print( "{0} = {1}".format( "outfile", codeml_params["outfile"] ), file=ctl )
	print( file=ctl )
	print( "{0} = {1}".format( "runmode", codeml_params["runmode"] ), file=ctl )
	print( "{0} = {1}".format( "seqtype", codeml_params["seqtype"] ), file=ctl )
	print( "{0} = {1}".format( "CodonFreq", codeml_params["CodonFreq"] ), file=ctl )
	print( "{0} = {1}".format( "model", codeml_params["model"] ), file=ctl )
	print( "{0} = {1}".format( "NSsites", codeml_params["NSsites"] ), file=ctl )
	print( "{0} = {1}".format( "fix_kappa", codeml_params["fix_kappa"] ), file=ctl )
	print( "{0} = {1}".format( "kappa", codeml_params["kappa"] ), file=ctl )
	print( "{0} = {1}".format( "fix_omega", codeml_params["fix_omega"] ), file=ctl )
	print( "{0} = {1}".format( "omega", codeml_params["omega"] ), file=ctl )
	print( "{0} = {1}".format( "ncatG", codeml_params["ncatG"] ), file=ctl )
print( "Finished Writing to: {0}".format( mini_ctl_path ) )
time.sleep( 1 )

with open( gene_dir_path + "omegas.lst", "w" ) as omfile:
	print( "\t".join( omega_lst[0] ), file=omfile )
	time.sleep(1)
	for i in tqdm(range( 0, aa_length-window_size+1 ), desc="Calculate omega", ascii=True):
	# for i in range( 0, aa_length-window_size+1 ):
		time.sleep(0.1) # sleep is redundant, but adding jut in case
		start = i
		stop = i+window_size
		temp_comparisons = []
		with open( temp_aln, "w" ) as outfile:
			for item in codon_dict.keys():
				print( ">{0}\n{1}".format( item, "".join(codon_dict[item][ start:stop ]) ), file=outfile )
				temp_comparisons.append( ">{0}\n{1}".format( item, "".join(codon_dict[item][ start:stop ]) ) )
		# os.system(  )

		run_line = "{0} {1} {2}".format( codeml_path, mini_ctl_path,  ">/dev/null 2>&1" )
		# print( "run_line: ", run_line )
		codeml_run = Popen([run_line], stdin=PIPE, shell=True)
		time.sleep(1)
		codeml_run.communicate( input = "\n".encode('utf-8') )
		# sys.exit( "in loop" )
		# if i == 10:
		# 	sys.exit( "test 10" )
		# os.system( "cat {0} | grep \" (dN/dS)\"".format( temp_mlc ) )
		# sys.exit( "test" )
		try:
			direct_output = subprocess.check_output("cat {0} | grep \"omega (dN/dS)\"".format( temp_mlc ), shell=True)
			omega = direct_output.strip().split()[-1].decode("utf-8")
			# print( "="*10,start,"="*10, stop,"="*10, omega  )
			# omega_lst.append( [start+1, stop+1, omega] )
			direct_output = subprocess.check_output("cat {0} | grep \"kappa (ts/tv)\"".format( temp_mlc ), shell=True)
			kappa = direct_output.strip().split()[-1].decode("utf-8")
			# print( "="*10,start,"="*10, stop,"="*10, kappa  )
			print( "\t".join([str(start+1), str(stop+1), str(omega), str(kappa), str(window_size)]), file=omfile )
			if float(omega) > 1:
				highly_evolvable[ "{0}:{1}".format( start, stop ) ] = temp_comparisons
			time.sleep(1)
			# print("here")
		except subprocess.CalledProcessError as e:
			raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
			sys.exit( "ERROR" )
			continue
		# sys.exit( "creation" )

# sys.exit("testing")


with open( "Tor_high_omegas.fna", "w" ) as outfile:
	for item in highly_evolvable.keys():
		print( ">REGION{0}\n{1}".format( item, highly_evolvable[item], file=outfile, end="\n\n" ) )
