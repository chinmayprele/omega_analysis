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
