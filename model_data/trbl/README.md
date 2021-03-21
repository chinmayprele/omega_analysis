## Primary Authors

- J. Yordy

## Gene Information

The gene _trbl_ functions in the Insulin Signaling Pathway's metabolism and cell growth through encoding the adapter protein "Tribbles" which negatively regulates signaling pathways in regards to the coordination of these biological processes.

The _Drosophila_ species annotated for the evolutionary rate analysis of _trbl_ were as follows:
	- _D. busckii_,
	- _D. arizonae_,
	- _D. persimilis_,
	- _D. suzukii_,
	- _D. willistoni_,
	- _D. hydei_,
	- _D. melanogaster_.

Then, Multiple Sequence Alignment using Clustal Omega within UGENE was done. dN/dS analysis was done using PAML's `codeml` on the whole gene, as well as sliding windows of size _n_ (dicated by the `omega_creator.py` script). dN/dS is relevant as it is the ratio of nonsynonymous substitutions to synonymous substitutions (nonsynonymous substitutions alter amino acids while synonymous substitutions do not alter amino acids). Foremost, this value is indicative of selective pressure, and it will be applied to analyze the evolutionary rate of trbl compared to the evolutionary rate of other genes in the Insulin Signaling Pathway with results being related to the conclusion of the Alvarez-Ponce paper "Network-level molecular evolutionary analysis of the insulin/TOR signal transduction pathway across 12 Drosophila genomes" in which genes more downstream exhibit more functional limitation.

Future directions of these analyses would be to investigate other genes along this pathway as well as the inclusion of more species to procure more precise data.

## Idiosyncrasies of Run

- Each windowed run should take about 5-8 seconds each with a few `time.sleep()` statements.
- However, a few runs ran for 130s+, which was unexpected.
	- Windows 6, 7, and 8 had this issue, but due to lack of constant monitoring, addition of other windows to this list is currently impossible.
	- There are certain runs that have very high omegas (and some with `NaN` omegas).
	- This might be the cause of the issue.
	- Furthermore, it could also be because "omegas (dN/dS)" could not be `grep`ped out of the file. This seems unlikely because near the `NaN`, the omegas are rather high, and the omegas could be artificially high. More analysis at this region needs to be incorporating more species.

There were also omega values that were rather high, and that brough up the omega value for the entire gene. These values were high due to indels present in the sequence, and the inability of `codeml` to distinguish between them, thus having infinite omega values (`nan`), and near-infinite omega values (`999.00`, which is the upper limit for `codeml`). To calculate the actual dN/dS of this gene, we had to omit the top 5 percentile of windowed values. "5" as a number was chosen randomly, but excision of these values landed us in an omega value we expected for this gene given its position in the pathway.

## References

- Alvarez-Ponce, David et al. “Network-level molecular evolutionary analysis of the insulin/TOR signal transduction pathway across 12 _Drosophila_ genomes.” Genome research vol. 19,2 (2009): 234-42. [doi:10.1101/gr.084038.108](https://genome.cshlp.org/content/19/2/234.short)
- Consortium, D. G., et al. (2007). "Evolution of genes and genomes on the Drosophila phylogeny." Nature 450(7167): 203-218.
- Okonechnikov K, Golosova O, Fursov M, the UGENE team. Unipro UGENE: a unified
	bioinformatics toolkit. Bioinformatics 2012 28: 1166-1167. [doi:10.1093/bioinformatics/bts091](https://pubmed.ncbi.nlm.nih.gov/22368248/)
- “Tribbles.” UniProt Consortium European Bioinformatics Institute Protein Information Resource SIB Swiss Institute of Bioinformatics, 2 Dec. 2020, [www.uniprot.org/uniprot/Q9V3Z1](www.uniprot.org/uniprot/Q9V3Z1).
