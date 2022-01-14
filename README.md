# A Phylodynamic Analysis of the Early COVID-19 Spread in the National Capital Region of the Philippines

Imari Joy C. Borda<sup>1</sup>, Sheryl Grace C. Buenaventura<sup>1</sup>,  Zython Paul T. Lachica<sup>1</sup>, Ivy Grace M. Panogalinog<sup>1,2</sup>, Ian Lorenzo Quibod<sup>3</sup>, Alexis Erich S. Almocera<sup>1,2</sup>, Lyre Anni E. Murao<sup>1,4</sup>, May Anne E. Mata<sup>1,2*</sup>, Ricardo C.H. del Rosario<sup>5</sup>

# Abstract

Computational analysis of genomic viral sequences has the potential to inform public health responses and to manage the spread of SARS-CoV-2. We sought to determine if phylodynamic modeling could describe the early local spread of the virus in the Philippines, and if the model could reflect the government’s control interventions. We applied previously published phylodynamic modeling tools to analyze 19 publicly available SARS-CoV-2 genomes sampled from 03 April to 18 July 2020 in the National Capital Region (NCR) of the Philippines. We estimate the date of the virus introduction, the effective reproduction number (R<sub>e</sub>), the basic reproduction number (R<sub>0</sub>), and the effective population size (N<sub>e</sub>). We found that the index case in NCR was on February 4, 2020 (95% Highest Posterior Density (HPD): November 18, 2019, to March 25, 2020) with R<sub>0</sub> of 2.20 (95% HPD: 1.74 to 2.80). Both the R<sub>e</sub> and the N<sub>e</sub> estimates suggest a progressive epidemic growth with an accelerated expansion on 20 May 2020 (95% HPD: 06 April 2020 – 08 June 2020), roughly the date of transition from  enhanced community quarantine (ECQ) to the more relaxed modified enhanced community quarantine (MECQ). Additionally, we estimate that there was underreporting in the total cases with only about 5% of the total cases being detected. Finally, the insights from different phylodynamic models coincide with each other. This demonstrates that even with limited genomic data, phylodynamic modeling can still provide useful insights about SARS-CoV-2, a virus  known to mutate rapidly whose evolution can be captured in just a short period. In summary, our results show a general view of the early stage spread of COVID-19 in NCR. Our study will serve as a conceptual framework for the use of SARS-CoV-2 genomic data in inferring epidemiological parameters to inform public health decisions in the Philippines.

# Content

* `data`: Contains the alignment of the sequences in FASTA format. The metadata was also included.
	* `ncr.fasta`: alignment used in the BEAST phylodynamic analysis.
	* `ncr-global.fasta`: alignment used in the creation of maximum likelihood tree.
	* `metadata-gisaid.tsv`: metadata of sequences downloaded from GISAID.
	* `metadata-ncbi.csv`: metadata of sequences downloaded fron NCBI.   
* `maximum-likelihood`: Contains the IQTree output of the maximum likelihood phylogeny.
	* `ncr-global.treefile`: ML tree which includes the NCR and global samples.
* `temporal-signal`: Contains the scripts and output files for calculating the temporal signal and plotting the corresponding root-to-tip divergence plot.
	* `ncr.treefile`: ML tree of the NCR samples using IQTree.
	* `genetic-distance.csv`: Datafile from TempEst which contains the genetic distance, dates, and residuals of the sequences.
	* `p-value.R`: Calculates the p-value of the root-to-stip regression. Also uses the function in `temsignalfunctions.R` and the R scripts inside the `functions` folder. These scripts were adapted from the work of Murray et al.
	* `p-value.RData`: Oytput from the `p-value.R`.
	* `regression-plot.R`: Plots the root-to-tip regression.
	* `temporal-signal.png`: Output from the `regression-plot.R`.
* `coalsky`: Contains BEAST xml file, R scripts, and outputs associated with the coalescent model. 
	* `ncr-coalsky.xml`: BEAST xml file for running the the coalescent skyline model.
	* `ncr-coalsky.log`: Output from the BEAST analysis.
	* `ncr-coalsky.csv`: Raw Ne estimates from the log file extracted using the Tracer tool.
	* `ncr-coalsky.R`: R script for plotting the Ne estimates.
	* `ncr-coalsky.png`: Output of the R script. This is the Ne plot. 
* `bdsky`: Contains BEAST xml file, R scripts, and outputs associated with the birth-deat skyline (BDSKY) model.
	* `ncr-bdsky-time-post.xml`: BEAST xml file for running the BDSKY model (posterior) for estimating the R<sub>e</sub> through time. The output of this is the `ncr-bdsky-time-post.log`. The corresponding prior run involves the xml file `ncr-bdsky-time-prior.xml` and its output is `ncr-bdsky-time-prior.log`.
	* `ncr-bdsky-time.R`: R script for plotting R<sub>e</sub> through time. The output is `ncr-bdsky-time.png`.
	* `ncr-bdsky-int-post.xml`: BEAST xml file for running the BDSKY model (posterior) for estimating the R<sub>e</sub> in relation to community quanrantine measures. The output of this is the `ncr-bdsky-int-post.log`. The corresponding prior runs involves the `ncr-bdsky-int-prior.xml` xml file and its output is `ncr-bdsky-int-prior.log`.
	* `ncr-bdsky-int.R`: R script for plotting R<sub>e</sub> in relation to community quanrantine measures. The output is `ncr-bdsky-int.png`.
	* `ncr-bdsky-maj-post.xml`: BEAST xml file for running the BDSKY model (posterior) for estimating the date of major change in R<sub>e</sub>. The output of this is the `ncr-bdsky-maj-post.log`. The corresponding prior runs involves the xml file `ncr-bdsky-maj-prior.xml` and its output is `ncr-bdsky-maj-prior.log`.
	* `ncr-bdsky-maj.R`: R script for plotting R<sub>e</sub> before and after the date of its major change. The output is `ncr-bdsky-maj.png`.
* `bdsir`: Contains BEAST xml file, R scripts, and outputs associated with the birth-death SIR (BDSIR) model. 
	* `ncr-bdsir-post.xml`: BEAST xml file for running the BDSIR model (posterior).
	* `ncr-bdsir-post.log`: Output from the BEAST analysis using `ncr-bdsir-post.xml`.
	* `ncr-bdsir-prior.xml`: BEAST xml file for running the BDSIR model (prior).
	* `ncr-bdsir-post.log`: Output from the BEAST analysis using `ncr-bdsir-prior.xml`.
	* `ncr-cases.csv`: Reported COVID-19 cases in the NCR, Philippines.
	* `ncr-bdsir-signal.R`: R script for plotting some of estimated parameters from the BDSIR model comparing the prior to the posterior.
	* `ncr-bdsir-signal.png`: Output from the `ncr-bdsir-signal.R`.
	* `loglist.txt`: List of log file (posterior only).
	* `ncr-bdsir-traj.R`: R script for plotting the SIR trajectories from the BDSIR model.
	* `ncr-bdsir-traj.png`: Output from the `ncr-bdsir-traj.R`.
