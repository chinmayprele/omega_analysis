## Gene Information

The gene _trbl_ functions in the Insulin Signaling Pathway's metabolism and cell growth through encoding the adapter protein "Tribbles" which negatively regulates signaling pathways in regards to the coordination of these biological processes.

The _Drosophila_ species annotated for the evolutionary rate analysis of trbl were as follows:
	- _D. busckii,_
	- _D. arizonae,_
	- _D. persimilis,_
	- _D. suzukii,_
	- _D. willistoni,_
	- _D. hydei,_
	- _D. melanogaster_.

Then, Multiple Sequence Alignment using Clustal Omega within UGENE was done. dN/dS analysis was done using PAML's `codeml` on the whole gene, as well as sliding windows of size _n_ (dicated by the `omega_creator.py` script). dN/dS is relevant as it is the ratio of nonsynonymous substitutions to synonymous substitutions (nonsynonymous substitutions alter amino acids while synonymous substitutions do not alter amino acids). Foremost, this value is indicative of selective pressure, and it will be applied to analyze the evolutionary rate of trbl compared to the evolutionary rate of other genes in the Insulin Signaling Pathway with results being related to the conclusion of the Alvarez-Ponce paper "Network-level molecular evolutionary analysis of the insulin/TOR signal transduction pathway across 12 Drosophila genomes" in which genes more downstream exhibit more functional limitation.

Future directions of these analyses would be to investigate other genes along this pathway as well as the inclusion of more species to procure more precise data.

## Primary Authors

- J. Yordy


### Idiosyncrasies of Run

- Each windowed run should take about 5-8 seconds each with a few `time.sleep()` statements.
- However, a few runs ran for 130s+, which was unexpected.
	- Windows 6, 7, and 8 had this issue, but due to lack of constant monitoring, addition of other windows to this list is currently impossible.
	- There are certain runs that have very high omegas (and some with `NaN` omegas).
	- This might be the cause of the issue.
	- Furthermore, it could also be because "omegas (dN/dS)" could not be `grep`ped out of the file. This seems unlikely because near the `NaN`, the omegas are rather high, and the omegas could be artificially high. More analysis at this region needs to be incorporating more species.
