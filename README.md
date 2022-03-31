# AncientEcoli
Scripts used for the paper: _A 16th Century Escherichia coli draft genome associated with an opportunistic bile infection_

The code and data starts after extracting the deduplicated mappings against the _E. coli_ genome. The steps used to get this point were the following:
	1. Trim and Merge the raw sequencing results using leeHom (with AdapterRemoval v2 identifying adapters)
	2. Pool lanes and digests
	3. Map against the human genome with BWA aln with the options `-n 0.01 -o 2 -l 16500`. Only reads with a minimum length of 30bp and minimum quality of 30 were kept
	4. Create an _E. coli_ pangenome using <Prokka> and <Roary>. The strains used to create the pan-genome are listed in Supplementary Data 1
	5. Reads that didn't map against the human genome are then mapped against this pan-genome. With the same settings as above.
	6. Reads that didn't successfully map against either the human or _E. coli_ genomes are finally mapped against _Klebsiella aerogenes_ (NC_015663.1).
The study can be found at: <DOI HERE>
The raw sequencing data can be found under the NCBI accession PRJNA810725
