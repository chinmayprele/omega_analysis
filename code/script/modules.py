"""
- Modules for omega_creator.py
"""

def codon_creator( aln_path, codons_path ):
	"""
	Creates a codons file from the alignment file.
	"""
	from Bio import SeqIO

	all_records = list(SeqIO.parse(aln_path, "fasta"))

	with open( codons_path, "w" ) as cfile:
		for record in all_records:
			sequence = str(record.seq)
			size_codon = 3
			codons = [sequence[i:i+size_codon] for i in range(0, len(sequence), size_codon)]
			print( "{0}\t{1}".format(
				record.id,
				" ".join(codons)
			), file=cfile )

	print( "Codons File Created!" )

	#
	# with open( aln_path ) as handle:
	#     for record in SeqIO.parse(handle, "clustal"):
	#         print(record.id)

def create_codon_file( codons_path, codon_dict ):
	"""
	creates a codon file from the ALN
	"""
	from tqdm import tqdm as tqdm

	with open( codons_path, "r") as infile:
		for line in tqdm(infile, desc="Reading in codons file", ascii=True):
			line = line.split()
			if line[0] not in codon_dict.keys():
				codon_dict[ line[0] ] = line[1:]
				aa_length = len( codon_dict[ line[0] ] )

def create_ctl( codeml_ctl_path, codeml_params ):
	"""
	creates codeml control file
	"""
	import time
	
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
