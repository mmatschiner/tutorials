# Divergence-Time Estimation with SNP Data

A tutorial on Bayesian divergence-time estimation with SNP data

## Summary

Besides [SVDQuartets](https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/) ([Chifman and Kubatko 2014](https://academic.oup.com/bioinformatics/article/30/23/3317/206559)), another method for phylogenetic inference with the multi-species-coalescent model based on SNP data is implemented in [SNAPP (SNP and AFLP Package for Phylogenetic analysis)](https://www.beast2.org/snapp/) ([Bryant et al. 2012](https://academic.oup.com/mbe/article/29/8/1917/1045283)), an add-on package for the program BEAST2. In principle, SNAPP is similar to the approach of StarBEAST2, only that each single SNP is considered as its own marker and gene trees are not separately inferred for each of these markers. Instead, SNAPP calculates the probability of the species tree without gene trees, by mathematically integrating over all possible gene trees. This approach reduces the parameter space of the model tremendously, which might be expected to reduce the time required for the analysis. Unfortunately, however, the mathematical integration over all possible gene trees is computationally very demanding; thus, SNAPP analyses are only feasible for a relatively small number of species.<br>Until recently, a limitation in the application of SNAPP has been that the branch lengths reported by SNAPP were in coalescent units rather than in units of time, and that these could not easily converted into time due to ascertainment bias in the SNP data. However, in the study of [Stange et al. (2018)](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy006/4827616), we tweaked the settings of SNAPP to include a strict-clock model that can be time calibrated based on the fossil record or on information from other phylogenies. In addition, we reduced the run time required for SNAPP analyses by linking all population sizes.


## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [SNP filtering](#filtering)
* [Divergence-time estimation with SNAPP](#snapp)


<a name="outline"></a>
## Outline

In this tutorial I am going to present how the BEAST2 add-on package SNAPP can be used for divergence-time estimation with SNP data. The settings of this model can not be specified through the BEAUti, but instead we will use the Ruby script [`snapp_prep.rb`](https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb) to generate input files for SNAPP. Differences between the species tree generated here based on SNP data and the species tree generated in tutorial [Bayesian Species-Tree Inference](../bayesian_species_tree_inference/README.md) will be discussed at the end of this tutorial.


<a name="dataset"></a>
## Dataset

The SNP data used in this tutorial is the filtered dataset used for species-tree inference with SVDQuartets in tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md). You can find more information about the origin of this dataset in the Dataset section of this other tutorial. In brief, the dataset has been filtered to include only bi-allelic SNPs with a low proportion of missing data, for the 26 samples of 13 cichlid species listed in the table below. Only SNPs mapping to chromosome 5 of the tilapia genome assembly ([Conte et al. 2017](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3723-5)) are included in the dataset, and these have been thinned so that no pair of SNPs is closer to each other than 100 bp.

<center>

| Sample ID | Species ID | Species name                  | Tribe         |
|-----------|------------|-------------------------------|---------------|
| IZA1      | astbur     | *Astatotilapia burtoni*       | Haplochromini |
| IZC5      | astbur     | *Astatotilapia burtoni*       | Haplochromini |
| AUE7      | altfas     | *Altolamprologus fasciatus*   | Lamprologini  |
| AXD5      | altfas     | *Altolamprologus fasciatus*   | Lamprologini  |
| JBD5      | telvit     | *Telmatochromis vittatus*     | Lamprologini  |
| JBD6      | telvit     | *Telmatochromis vittatus*     | Lamprologini  |
| JUH9      | neobri     | *Neolamprologus brichardi*    | Lamprologini  |
| JUI1      | neobri     | *Neolamprologus brichardi*    | Lamprologini  |
| KHA7      | neochi     | *Neolamprologus chitamwebwai* | Lamprologini  |
| KHA9      | neochi     | *Neolamprologus chitamwebwai* | Lamprologini  |
| IVE8      | neocra     | *Neolamprologus crassus*      | Lamprologini  |
| IVF1      | neocra     | *Neolamprologus crassus*      | Lamprologini  |
| JWH1      | neogra     | *Neolamprologus gracilis*     | Lamprologini  |
| JWH2      | neogra     | *Neolamprologus gracilis*     | Lamprologini  |
| JWG8      | neohel     | *Neolamprologus helianthus*   | Lamprologini  |
| JWG9      | neohel     | *Neolamprologus helianthus*   | Lamprologini  |
| JWH3      | neomar     | *Neolamprologus marunguensis* | Lamprologini  |
| JWH4      | neomar     | *Neolamprologus marunguensis* | Lamprologini  |
| JWH5      | neooli     | *Neolamprologus olivaceous*   | Lamprologini  |
| JWH6      | neooli     | *Neolamprologus olivaceous*   | Lamprologini  |
| ISA6      | neopul     | *Neolamprologus pulcher*      | Lamprologini  |
| ISB3      | neopul     | *Neolamprologus pulcher*      | Lamprologini  |
| ISA8      | neosav     | *Neolamprologus savoryi*      | Lamprologini  |
| IYA4      | neosav     | *Neolamprologus savoryi*      | Lamprologini  |
| KFD2      | neowal     | *Neolamprologus walteri*      | Lamprologini  |
| KFD4      | neowal     | *Neolamprologus walteri*      | Lamprologini  |

</center>


<a name="requirements"></a>
## Requirements

* **BEAST2:** If you followed other tutorials in this collection, you likely have the BEAST2 package, including BEAST2 itself and TreeAnnotator, installed already. If not, you can downloaded these from the BEAST2 website [https://www.beast2.org](https://www.beast2.org). As all these programs are written in Java, compilation is not required, and all programs should work on Mac OS X, Linux, and Windows. I recommend downloading the program versions that include the Java Runtime Environment, which may prevent conflicts with Java versions that may already be installed on your machine.<br>

* **Tracer:** Like BEAST2, you probably have Tracer installed already if you did other tutorials of this collection. If not, you can download Tracer for Mac OS X, Linux, or Windows from [http://tree.bio.ed.ac.uk/software/tracer/](http://tree.bio.ed.ac.uk/software/tracer/). The file with the extension `.dmg` is for Mac OS X, the one with the extension `.tgz` is for Linux, and the Windows version is the file ending in `.zip`.
		
* **bcftools:** If you did tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md), you probably have [bcftools](http://www.htslib.org/doc/bcftools.html) ([Li 2011](https://academic.oup.com/bioinformatics/article/27/21/2987/217423)) installed already. If not, you can find downloads and installation instructions for Mac OS X and Linux at the [HTSlib download webpage](http://www.htslib.org/download/). Installation on Windows is apparently not possible. If you should fail to install bcftools, you could skip the optional [SNP filtering](#filtering) steps in the first part of the tutorial and still run the [Divergence-time estimation with SNAPP](#snapp).



<a name="filtering"></a>
## SNP filtering

Based on simulations, we tested the performance of SNAPP with a range of datasets in [Stange et al. (2018)](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy006/4827616). We found that the run time of SNAPP increases more or less linearly with the number of SNPs included in the dataset, but that that number of samples used per species has an even stronger influence on run time. However, we also found that accurate and strongly supported species can be obtained with dataset containing only around 1,000 SNPs for just a single sample per species. Thus, before preparing the input file for SNAPP, we will first reduced the dataset that was already filtered in tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md). This further filtering will again be done in bcftools, but if you failed to install bcftools, you could skip the following steps and continue below with [Divergence-time estimation with SNAPP](#snapp)

* Copy file [`NC_031969.f5.sub4.vcf`](data/NC_031969.f5.sub4.vcf) to your current analysis directory, either from the directory that you used for tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md) (if you ran that other tutorial) or by downloading this file with the link.

* As mentioned above and in tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md), the dataset in file `NC_031969.f5.sub4.vcf` contains SNP data for two samples per species. To reduce the run time for SNAPP, we are going to generate another version of the same dataset that is reduced to only a single sample per species. To exclude the samples "IZA1", "AXD5", "JBD5", "JUI1", "KHA9", "IVF1", "JWH1", "JWG8", "JWH3", "JWH5", "ISA6", "IYA4", and "KFD4" (these were selected because they have more missing data than the other samples of each species) from a new file named `NC_031969.f5.sub5.vcf`, use the following command:

		bcftools view -s ^IZA1,AXD5,JBD5,JUI1,KHA9,IVF1,JWH1,JWG8,JWH3,JWH5,ISA6,IYA4,KFD4 -o NC_031969.f5.sub5.vcf NC_031969.f5.sub4.vcf
		
* The reduction of samples in the previous step might have led to some sites becoming monomorphic for either the reference or the alternate allele. Some sites might also have no data left for one or more of the species, given that these are now represented by only one instead of two samples. To exclude these sites again with bcftools, use the following command:

		bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.0' -o NC_031969.f5.sub6.vcf NC_031969.f5.sub5.vcf
	
	**Question 1:** How many SNPs are now left in file `NC_031969.f5.sub6.vcf`? [(see answer)](#q1)


<a name="snapp"></a>
## Divergence-time estimation with SNAPP

In this part of the tutorial, we are going to prepare the XML input file for SNAPP with the Ruby script `snapp_prep.rb` so that the model settings of [Stange et al. (2018)](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy006/4827616) are used for divergence-time estimation with SNAPP. In this model, the population sizes of all species are linked and the substitution rate (which is assumed to be according to a strict clock) can be calibrated with age constraints on one or more divergence events.

* First, download the script [`snapp_prep.rb`](https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb) from github, either by clicking on the link or by using the following command:

		wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb
		
* Then, to see the options available for generating SNAPP input files with `snapp_prep.rb`, have a look at the help text of the script by using the following command:

		ruby snapp_prep.rb -h
		
	You'll see that `snapp_prep.rb` accepts input either in Phylip (with option "-p") or VCF (option "-v") format. In addition, the script requires a table file (option "-t") in which samples are assigned to species, and a constraint file (option "-c") in which age constraints for selected divergence events are specified. None of the other options is required, but these can be useful if you need to specify a starting tree (option "-t") because SNAPP fails to find a suitable starting tree itself, if you want to increase or decrease the default length (500,000 iterations) of the MCMC chain (option "-l"), or if you want to limit the dataset to a certain number of SNPs to reduce the run time (option "-m").
	
* Examples for the table and constraint files can be found on the [github repository for `snapp_prep.rb`](https://github.com/mmatschiner/snapp_prep). Have a look at these example files; the example table file is named [`example.spc.txt`](https://github.com/mmatschiner/snapp_prep/blob/master/example.spc.txt) and the example constraint file is named [`example.con.txt`](https://github.com/mmatschiner/snapp_prep/blob/master/example.con.txt).

* To prepare the table file assigning samples to species, write the following text to a new file named `specimens.txt`:

		species	specimens
		astbur	IZC5
		altfas	AUE7
		telvit	JBD6
		neobri	JUH9
		neochi	KHA7
		neocra	IVE8
		neogra	JWH2
		neohel	JWG9
		neomar	JWH4
		neooli	JWH6
		neopul	ISB3
		neosav	ISA8
		neowal	KFD2

* For the constraint file, we'll need to specify at least one age constraint for a divergence among the 13 cichlid species. For this, we can refer to the results of the analysis with the multi-species coalescent model in tutorial [Bayesian Species-Tree Inference](../bayesian_species_tree_inference/README.md). In that tutorial, the divergence of *Astatotilapia burtoni* ("astbur") and the *Neolamprologus* species was estimated at 6.2666 Ma with a 95% HPD interval from 4.8225 to 7.9049 Ma (if you'ld like to look up these results you can find them in file [`starbeast_species.tre`](data/starbeast_species.tre)). We can approximate this mean and confidence interval with a lognormal distribution centered at 6.2666 Ma that has a standard deviation (in real space) of 0.13. In addition, we could also constrain the first divergence among the species *Neolamprologus brichardi* ("neobri"), *Neolamprologus gracilis* ("neogra"), *Neolamprologus marunguensis* ("neomar"), and *Neolamprologus olivaceus* ("neooli") according to the results of this earlier analysis; however, it is likely that this divergence event would be more precisely estimated with the current SNP dataset. Thus, we are only going to constrain the age of a single node: the divergence of *Astatotilapia burtoni* ("astbur") from the species of the tribe Lamprologini.<br>To write this single age constraint to a constraint file, write the following text to a new file named `constraints.txt`:

		lognormal(0,6.2666,0.13)  crown astbur,altfas,telvit,neobri,neochi,neocra,neogra,neohel,neomar,neooli,neopul,neosav,neowal
		
	Note that the above line contains three character strings that are delimited with tab symbols. The first of the three character strings specifies that a lognormal distribution with an offset of 0, a mean of 6.2666, and a standard deviation (in real space) of 0.13 should be used for the constraint. The second character string ("crown") specifies that the crown divergence of the clade should be constrained rather than the stem divergence (this would be the time at which the clade diverged from its sister group). Finally, the third character string simply lists all species included in the constrained clade. Because the constraint used here applies to the very first divergence of the phylogeny (the root), all species names are listed in the third character string. **Question 2:** Can you think of a way how the exact same age calibration could be specified differently? [(see answer)](#q2)

* With the table file and the constraints file prepared, we can now prepare the input file for SNAPP with the script `snapp_prep.rb`. To limit the dataset to 1,000 randomly selected SNPs and to set a chain length of 100,000 MCMC iterations, we'll use the options "-m 1000" and "-l 100000", respectively. Thus, use the following command to generate the XML input file for SNAPP:

		ruby snapp_prep.rb -v NC_031969.f5.sub6.vcf -t specimens.txt -c constraints.txt -m 1000 -l 100000
		
	You may notice that the chain length of 100,000 MCMC iterations is extremely short compared to those used in the BEAST2 analyses of other tutorials. Using much shorter chain lengths with SNAPP than BEAST2 is quite common, given that the SNAPP model has far fewer model parameters than most models used in other BEAST2 analyses, and that each individual iteration is much slower with SNAPP due to the high computational demand of the integration over all possible gene trees at each SNP.
		
* The script `snapp_prep.rb` should have written a file named `snapp.xml`. You could open that file in a text editor and read some of the annotations if you'ld like to know more about the settings that we are about to use with SNAPP.
		
* Before we can run SNAPP, we still need to install this add-on package for BEAST2. As with BEAST2 add-on packages used in other tutorials, this can be done with the program BEAUti. Open BEAUti, then click on "Manage Packages" in BEAUti's "File" menu. This will open a new window for the BEAST2 Package Manager. In that window, select "SNAPP" and click "Install/Upgrade" as shown in the screenshot below. Then, close BEAUti again.<p align="center"><img src="img/beauti1.png" alt="BEAUti" width="700"></p>

* To "run SNAPP", we actually run BEAST2, either using the GUI or the command-line version. In both cases, we can speed up the analysis by using multiple CPUs, as SNAPP analyses are highly parallelizable. Thus, if you have four CPUs available on your machine and use all of them for the SNAPP analysis, this analysis should take only about a fourth of the time that would be required with a single CPU. If you should run SNAPP analyses on a server, you might even be able to use tens of CPUs simultaneously, which would shorten SNAPP's run times tremendously. To start the analysis for example on four CPUs with the GUI version of BEAST2, select "4" from the drop-down menu next to "Thread pool size" as shown in the next screenshot.<p align="center"><img src="img/beast1.png" alt="BEAST" width="500"></p>If you would use the command-line version of BEAST instead, you could start the analysis with four CPUs using the following command:

		/Applications/Beast/2.5.0/bin/beast -threads 4 snapp.xml  

	(if BEAST2 installation is not in `/Applications/Beast/2.5.0` on your machine, you'll have to replace this part of the command with the actual path to your BEAST2 installation).
	
	Depending on the number of CPUs that you have available for SNAPP this analysis will probably take between 2-3 hours (with four CPUs) and 10-12 hours (with a single CPU). If you do not want to wait this long before continuing the tutorial, you can find the results of my analysis in files [`snapp.log`](res/snapp.log) and [`snapp.trees`](res/snapp.trees).
	
* Once the SNAPP analysis has completed or if you decided to use the results of my analysis instead of generating your own, open file `snapp.log` in Tracer. The Tracer window should then display run statistics similar (or identical if you used the results of my analysis) to those shown in the next screenshot.<p align="center"><img src="img/tracer1.png" alt="Tracer" width="700"></p>You'll see that some low ESS values indicate that the MCMC chain has not yet reached stationarity and that the analysis should ideally have been performed with more MCMC iterations. For our interpretation here, however, we'll assume that the degree of stationarity is sufficient. What you also should notice is that the list of parameters on the left-hand side of the window is now much shorter than it usually is with results of BEAST2 analyses. In fact, only three parameters are shown: The speciation rate ("lambda"), the age of the root of the species tree ("treeHeightLogger"), and the substitution rate ("clockRate"). Note that the substitution rate is not comparable to a genome-wide rate due to ascertainment bias in the SNP dataset, even though SNAPP by default applies an ascertainment-bias correction (see [Bryant et al. 2012](https://academic.oup.com/mbe/article/29/8/1917/1045283) for details). However, one model parameter is not included in the log output yet, namely the population size. Recall that by using script `snapp_prep.rb` to generate the XML file, we implemented a model in which the population sizes of all branches are set to be identical to each other. Thus, if the population sizes had been included in the log file, this file would contain a large number of columns with identical information. To avoid this, the output of the population sizes to the log file has been disabled by `snapp_prep.rb`. However, the population size estimates are still available because they were instead written to the tree file `snapp.trees`, and we can now add them to the log file `snapp.log` using the Ruby script `add_theta_to_log.rb`. Thus, download the script [`add_theta_to_log.rb`](https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/add_theta_to_log.rb) from the github repository for `snapp_prep.rb`, either by clicking the link or with the following command:

		wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/add_theta_to_log.rb
		
* Have a look at the available options for `add_theta_to_log.rb`:

		ruby add_theta_to_log.rb -h
		
	You'll see that this script requires four command-line arguments: The names of the log and tree input files, the name of an output file, as well as an estimate for the generation size. The latter is required to calculate an estimate of the effective population size, for which `add_theta_to_log.rb` uses the equation <i>N</i><sub>e</sub> = Theta &div; (4 &times; <i>r</i> &div; <i>n</i><sub>g</sub>), where <i>r</i> is the substitution-rate estimate (= the rate of the strict clock) and <i>n</i><sub>g</sub> is the number of generations per time unit (more details are given in [Stange et al. 2018](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy006/4827616)).
	
* We'll here assume (as in other tutorials) that the generation time of cichlids is three years, and we'll name the output file `snapp_w_popsize.log`. Thus use the following command to run script `add_theta_to_log.rb`:

		ruby add_theta_to_log.rb -l snapp.log -t snapp.trees -g 3 -o snapp_w_popsize.log

* Then, open file `snapp_w_popsize.log` again in Tracer. You should see that "theta" and "population_size" have now been added to the list of parameters, as shown in the next screenshot.<p align="center"><img src="img/tracer2.png" alt="Tracer" width="700"></p>

* Select both "clockRate" and "theta" from the list of parameters and click the tab for "Joint-Marginal" at the top right. You should then see that both of these parameters are highly correlated, as shown in the next screenshot.<p align="center"><img src="img/tracer3.png" alt="Tracer" width="700"></p>Note that both of these estimates should not be taken as being representative for the whole genome due to ascertainment bias. Because invariant sites are not included in a SNP dataset, this SNP dataset necessarily appears more variable than the whole-genome dataset from which it was derived, and thus the substitution rate of the SNP data appears to be higher than the rate would be if it was averaged over the entire genome. Nevertheless, our simulations in [Stange et al. (2018)](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy006/4827616) have shown that the population-size estimates of SNAPP are actually reliable.

* Select "population_size" from the list of parameters and click on the "Estimate" tab to see summary information and a histogram for the population-size estimate, as shown in the next screenshot.<p align="center"><img src="img/tracer4.png" alt="Tracer" width="700"></p> **Question 3 (if you already did tutorial Bayesian Species-Tree Inference):** How does this estimate of the population size of Lake Tanganyika cichlid fishes compare to that assumed in tutorial [Bayesian Species-Tree Inference](../bayesian_species_tree_inference/README.md)? [(see answer)](#q3)

* Next, open the file `snapp.trees` in the software Densitree from the BEAST2 package to visualize the full set of posterior trees sampled by SNAPP, as shown in the next screenshot.<p align="center"><img src="img/densitree1.png" alt="DensiTree" width="600"></p> You should see that not all posterior trees share the same topology, indicating remaining uncertainty in the relationships of the 13 cichlid species. In particular, the relationships of *Telmatochromis vittatus* ("telvit") appear ambiguous, as this species is placed next to *Neolamprologus walteri* ("neowal") and *Neolamprologus chitamwebwai* ("neochi") in some of the posterior trees, but apparently closer to the remaining *Neolamprologus* species or ancestral to all of them in the other posterior trees. The relationships among the five species *Neolamprologus brichardi* ("neobri"), *Neolamprologus olivaceous* ("neooli"), *Neolamprologus pulcher* ("neopul"), *Neolamprologus helianthus* ("neohel"), and *Neolamprologus gracilis* ("neogra") also appear uncertain.

* Also quantify the posterior probabilities of clades as node support in a maximum-clade-credibility tree using TreeAnnotator. To do so, open `snapp.trees` in TreeAnnotator, set the burnin percentage to 10, choose "Mean heights" from the drop-down menu next to "Node heights", select "snapp.trees" as the input file, and name the output file "snapp.tre", as shown in the next screenshot. Then, click "Run".<p align="center"><img src="img/treeannotator1.png" alt="TreeAnnotator" width="500"></p>

* Finally, open file [`snapp.tre`](res/snapp.tre) in FigTree and display the "posterior" node support values as node labels, as shown in the screenshot below.<p align="center"><img src="img/figtree1.png" alt="FigTree" width="600"></p> The posterior probabilities for the different clades support the interpretation based on the Densitree plot made above: The position of *Telmatochromis vittatus* ("telvit") is uncertain, as are the relationships among the five species *Neolamprologus brichardi* ("neobri"), *Neolamprologus olivaceous* ("neooli"), *Neolamprologus pulcher* ("neopul"), *Neolamprologus helianthus* ("neohel"), and *Neolamprologus gracilis* ("neogra"). However, the visualization in FigTree also shows that a clade comprising the three species *Neolamprologus brichardi* ("neobri"), *Neolamprologus olivaceous* ("neooli"), and *Neolamprologus pulcher* ("neopul") is actually well supported; this has not been apparent in DensiTree. **Question 4 (if you already did tutorial Species-Tree Inference with SNP Data):** Is the species-tree topology estimated with SNAPP concordant with that estimated with SVDQuartets in tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md)? [(see answer)](#q4) **Question 5 (if you already did tutorial Bayesian Species-Tree Inference):** How does the divergence-time estimate for *Neolamprologus marunguensis* ("neomar"), *Neolamprologus gracilis* ("neogra"), *Neolamprologus brichardi* ("neobri"), and *Neolamprologus olivaceous* ("neooli") compare to that obtained with StarBEAST2 in tutorial [Bayesian Species-Tree Inference](bayesian_species_tree_inference/README.md)?

Despite the remaining uncertainty in the relationships among the *Neolamprologus* species, the SNAPP analysis has allowed us to estimate an average population size for species of the tribe Lamprologini, and it has improved the estimates of divergence times within the group. The population-size estimate will be useful in the inference of introgression events in tutorial [Bayesian Analysis of Reticulate Evolution](../bayesian_analysis_of_reticulate_evolution/README.md), and the improved divergence times will be used for time calibration in tutorial [Analysis of Introgression with Chromosome-Length Alignments](../analysis_of_introgression_with_chromosome_length_alignments/README.md).

<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** The reduction from two samples to one sample per species led to a large decrease in the number of bi-allelic SNPs. Only 18,195 SNPs should now remain in the filtered dataset of file `NC_031969.f5.sub6.vcf`, which nevertheless is more than sufficient for a phylogenetic analysis with SNAPP. To see the number of SNPs in this VCF file, you can use the the following command:

		bcftools view -H NC_031969.f5.sub6.vcf | wc -l
		
		
<a name="q2"></a>

* **Question 2:** Because the crown divergence of a clade is identical to the stem divergences of the two descending clades, we could also specify the same age constraint with the following line:

		lognormal(0,6.2666,0.13)  stem altfas,telvit,neobri,neochi,neocra,neogra,neohel,neomar,neooli,neopul,neosav,neowal
		
	Note that "stem" is now specified as the second character string instead of "crown" and that "astbur" is now missing from the list of species in the third character string.
	

<a name="q3"></a>

* **Question 3:** Recall that in tutorial [Bayesian Species-Tree Inference](../bayesian_species_tree_inference/README.md), the population-size parameter was fixed in the analysis with the multi-species coalescent model according to an estimate published by [Meyer et al. (2017)](https://academic.oup.com/sysbio/article/66/4/531/2670093). In this study, [Meyer et al. (2017)](https://academic.oup.com/sysbio/article/66/4/531/2670093) reported that "effective population sizes (<i>N</i><sub>e</sub>) estimated with the multispecies coalescent model ranged between 3.6 &times; 10<sup>4</sup> and 8.1 &times; 10<sup>5</sup>, assuming a mean generation times of 3 years for cichlid fishes". Most of these population-size estimates in [Meyer et al. (2017)](https://academic.oup.com/sysbio/article/66/4/531/2670093) were around 3.3 &times; 10<sup>5</sup>, therefore this value was used in tutorial [Bayesian Species-Tree Inference](../bayesian_species_tree_inference/README.md). In comparison to this population size assumed in the other tutorial, the current population-size estimate based on the SNAPP analysis, with a mean of 9.3 &times; 10<sup>4</sup> and a 95% HPD interval ranging from 6.5 &times; 10<sup>4</sup> to 1.2 &times; 10<sup>5</sup>, is much smaller.


<a name="q4"></a>

* **Question 4:** The species trees estimated with SNAPP and SVDQuartets are not fully concordant. If you'ld like to look up the result of the SVDQuartets analysis, open file [`NC_031969.f5.sub4.parts.tre`](data/NC_031969.f5.sub4.parts.tre) in FigTree, and display "label" (these are the bootstrap support values) as node labels. In contrast to SNAPP, SVDQuartets had strongly supported a position of *Telmatochromis vittatus* ("telvit") outside of a monophyletic clade comprising all representatives of the genus *Neolamprologus*. In addition, SVDQuartets had placed *Neolamprologus helianthus* ("neohel") as the sister of species trio *Neolamprologus brichardi* ("neobri"), *Neolamprologus olivaceous* ("neooli"), and *Neolamprologus pulcher* ("neopul") with similarly strong support, whereas *Neolamprologus helianthus* ("neohel") appeared as the sister of *Neolamprologus gracilis* ("neogra") in the SNAPP analysis. Given that in both cases of conflict between the two analyses, the SVDQuartets analyses had inferred strong support based on bootstrapping whereas SNAPP only assigned moderate posterior probabilities, it is likely that the relationships estimated by SVDQuartets is the correct one. The lower support values in the SNAPP analysis are not surprising, though, given that only a small fraction of the data used with SVDQuartets (over 63,000 SNPs) has been used in the analysis with SNAPP (1,000 SNPs). It may be expected that using a larger number of SNPs with SNAPP would increase the reliability of the resulting species tree.


<a name="q5"></a>

* **Question 5:** Recall that the divergence time of the four species *Neolamprologus marunguensis* ("neomar"), *Neolamprologus gracilis* ("neogra"), *Neolamprologus brichardi* ("neobri"), and *Neolamprologus olivaceous* ("neooli") was estimated around 0.93 Ma in the analysis with StarBEAST2 in tutorial [Bayesian Species-Tree Inference](bayesian_species_tree_inference/README.md). To look up these results, you could open file [`starbeast_species.tre`](data/starbeast_species.tre) in FigTree and display "Node ages" as node labels. In contrast, the first divergence of the clade comprising these four species was estimated with SNAPP at a more recent time, around 0.60 Ma, as shown in the screenshot below. Moreover, the 95% HPD interval for this age estimate ranges from 0.41 to 0.78 Ma (select "height\_95%\_HPD" instead of "Node ages" from the drop-down menu to see these) and therefore does not even include the older age estimate of the StarBEAST2 analysis. In contrast the 95% HPD interval for this divergence time in the StarBEAST2 analysis was wider and ranged from 0.42 to 1.44 Ma, which thus includes almost the entire 95% HPD interval of the SNAPP analysis. This indicates that the two analyses do not disagree with each other regarding this divergence time but only that the StarBEAST2 analysis was not able to estimate it as precisely as the SNAPP analysis.<p align="center"><img src="img/figtree2.png" alt="FigTree" width="600"></p>