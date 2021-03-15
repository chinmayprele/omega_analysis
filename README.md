# omega_analysis

We will be studying the dN/dS ratio for sliding windows in multiple genes within the _Drosophila_ phylogeny.

We will be using `codeml` from the [Phylogenetic Analysis by Maximum Likelihood (PAML)](http://abacus.gene.ucl.ac.uk/software/paml.html) in order to do this.

CPR Local Path: `/Users/rele.c/Downloads/omega_analysis`

### To Add Files
- Within `model_data`, create a directory with the name of your gene, make sure it is unique.
- Make the following subdirectories within your gene:
	- `ctl_files`: empty, needed for run
	- `images`: to contain TIFFs of your images
	- `raw_data`: to contain raw data of your models.
- The `raw_data` for each will contain two file, the MSA `aln` from [UGENE](http://ugene.net), and a subset of the tree as a `nwk` from [iTOL](https://itol.embl.de).
	- Please name your `aln` file as: `<gene_name>.msa.aln`. _Within the file, make sure to remove the name of the gene and the transcript from the header_. The only name that should be there is the name of the assembly.
 	- Please name your `nwk` file as: `<gene_name>.nwk`. We might have to change the format  of the `nwk`, but can do that if we need to.
	If these files are named wrong, the code __will not run__.
- Please add a `README.md` as well. It should contain your names, list of all species analyzed, as well as some information for each gene.

### To Run

- If you plan to run this, you will have to change the paths of the files in `omega_analysis.py` to match your local directory.
	- [ ] I may make them relative, but for now, they are absolute.
- Go to the directory with `omega_analysis.py`, and run it using `python3 omega_analysis.py`.
- Then select the gene from the list that pops up. Hit "Enter".
- You may have to wait until the progress bar for the window analysis pops up, which should take between 2-5 minutes.

#### Changing Parameters
- If needed, you can change the size of the sliding windows in that file before the run.
- You may also change the `codeml.ctl` parameter options in the same file. Please refer to the [PAML Documentation](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf) for more information.

#### After the Run
- You will have to clean up unnecessary files, specifically:
	- From `code/script`
		- `2AA.t`
		- `2ML.dN`
		- `2ML.dS`
		- `2ML.t`
		- `2NG.dS`
		- `2NG.t`
		- `4fold.nuc`
		- `Inf `
		- `rst`
		- `rst1`
		- `rub`
	- From `model_data/<gene>`
		- `<gene>.mlc`
		- `temp.aln`
		- `temp.mlc`
- The script should still work without the cleanup, since it is overriding the data for each run, but it is good to remove it regardless.

### System Config
Running on:
- MacBook Pro (13-inch, 2019, Four Thunderbolt 3 ports)
- _Processor_: 2.4 GHz Quad-Core Intel Core i5
- _Memory_: 16 GB 2133 MHz LPDDR3
- _Graphics_: Intel Iris Plus Graphics 655 1536 MB
- _Shell_: `bash`
- _OS_: macOS Big Sur (Version 11.2.2)

### Acknowledgements

We would like to thank [D'Andrew L. Harrington](https://github.com/ContinuumDLH) who helped us set parameters and get `codeml` running.
