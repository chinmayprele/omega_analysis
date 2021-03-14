# omega_analysis

We will be studying the dN/dS ratio for sliding windows in multiple genes within the _Drosophila_ phylogeny.

We will be using `codeml` from the [Phylogenetic Analysis by Maximum Likelihood (PAML)](http://abacus.gene.ucl.ac.uk/software/paml.html) in order to do this.

CPR Local Path: `/Users/rele.c/Box/00_Personal/reed_lab/URCA/2021_spring/omega_analysis`

### To add files:
- Within `model_data`, create a directory with the name of your gene, make sure it is unique.
- Make the following subdirectories within your gene:
	- `ctl_files`: empty, needed for run
	- `images`: to contain TIFFs of your images
	- `raw_data`: to contain raw data of your models.
- The `raw_data` for each will contain two file, the MSA `aln` from [UGENE](http://ugene.net), and a subset of the tree as a `nwk` from [iTOL](https://itol.embl.de).
	- Please name your `aln` file as: `<gene_name>.msa.aln`. Within the file, make sure to remove the name of the gene and the transcript from the header. The only name that should be ther eis the name of the assembly.
 	- Please name your `nwk` file as: `<gene_name>.nwk`. We might have to change the format  of the `nwk`, but can do that if we need to.
	If these files are named wrong, the code __will not run__.
- Please add a `README.md` as well. It should contain your names, list of all species analyzed, as well as some information for each gene.

### System Config
Running on:
- MacBook Pro (13-inch, 2019, Four Thunderbolt 3 ports)
- __Processor__: 2.4 GHz Quad-Core Intel Core i5
- __Memory__: 16 GB 2133 MHz LPDDR3
- __Graphics__: Intel Iris Plus Graphics 655 1536 MB

### Acknowledgements

We would like to thank [D'Andrew L. Harrington](https://github.com/ContinuumDLH) who helped us set parameters and get `codeml` running.
