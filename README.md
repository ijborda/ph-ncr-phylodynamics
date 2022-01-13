# A Phylodynamic Analysis of the Early COVID-19 Spread in the National Capital Region of the Philippines

# Organization

* `data`: Contains the alignment of the sequences in FASTA format. The metadata was also included.
	* `ncr.fasta`: alignment used in the BEAST phylodynamic analysis.
	* `ncr-global.fasta`: alignment used in the creation of maximum likelihood tree.
	* `metadata-gisaid.tsv`: metadata of sequences downloaded from GISAID.
	* `metadata-ncbi.csv`: metadata of sequences downloaded fron NCBI.   
* `maximum-likelihood`: Contains the IQTree output of the maximum likelihood phylogeny.
	* `ncr-global.treefile`: ML tree which includes the NCR and global samples.
* `temporal-signal`: Contains the scripts and output files for calculating the temporal signal and plotting the corresponding root-to-tip divergence plot.
	* `ncr.treefile`: ML tree of the NCR samples using IQTree.
	* `genetic-distance`: Datafile from TempEst which contains the genetic distance, dates, and residuals of the sequences.
	* `p-value.R`: Calculates the p-value of the root-to-stip regression. Also uses the function in `temsignalfunctions.R` and the R scripts inside the `functions` folder. These scripts were adapted from the work of Murray et al.
	* `p-value.RData`: Oytput from the `p-value.R`.
	* `regression-plot.R`: Plots the root-to-tip regression.
	* `temporal-signal.png`: Output from the `regression-plot.R`.
