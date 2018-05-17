# Bayesian Species-Tree Inference

A tutorial on Bayesian inference of time-calibrated species trees

## Summary

XXX

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Bayesian species-tree inference with StarBEAST2](#starbeast2)
* [Bayesian species-tree inference with concatenation](#concatenation)
* [Comparing divergence times estimated with StarBEAST2 and concatenation](#comparison)


<a name="outline"></a>
## Outline

In this tutorial, XXX


<a name="dataset"></a>
## Dataset

As in tutorial [Maximum-Likelihood Species-Tree Inference](../ml_species_tree_inference/README.md), the dataset used here is the set of alignments for 72 genes produced in tutorial [Ortholog Detection](../ortholog_detection/README.md). This dataset includes sequences for eleven cichlid species, two of which represent Neotropical cichlids while the remaining nine species are from Africa. The focus of the taxon set is on cichlids of the rapid radiations in the East African lakes Tanganyika, Malawi, and Victoria. The table below lists all species included in the set of alignments. Note that the sequence data for *Ophthalmotilapia ventralis* were extracted from a transcriptome assembly whereas genome assemblies were available for all other species.

<center>

| ID     | Species                         | Tribe          | Distribution    |
|--------|---------------------------------|----------------|-----------------|
| ampcit | *Amphilophus citrinellus*       | Heroini        | Neotropics      |
| andcoe | *Andinoacara coeruleopunctatus* | Cichlasomatini | Neotropics      |
| orenil | *Oreochromis nilotiucs*         | Oreochromini   | African rivers  |
| ophven | *Ophthalmotilapia ventralis*    | Ectodini       | Lake Tanganyika |
| astbur | *Astatotilapia burtoni*         | Haplochromini  | Lake Tanganyika |
| metzeb | *Metriaclima zebra*             | Haplochromini  | Lake Malawi     |
| punnye | *Pundamilia nyererei*           | Haplochromini  | Lake Victoria   |
| neobri | *Neolamprologus brichardi*      | Lamprologini   | Lake Tanganyika |
| neomar | *Neolamprologus marunguensis*   | Lamprologini   | Lake Tanganyika |
| neogra | *Neolamprologus gracilis*       | Lamprologini   | Lake Tanganyika |
| neooli | *Neolamprologus olivaceous*     | Lamprologini   | Lake Tanganyika |

</center>

<a name="requirements"></a>
## Requirements

* **BEAST2:** The BEAST2 package, including BEAUti, BEAST2 itself, TreeAnnotator, and other tools can be downloaded from the BEAST2 website [https://www.beast2.org](https://www.beast2.org). As all these programs are written in Java, compilation is not required, and all programs should work on Mac OS X, Linux, and Windows. I recommend downloading the program versions that include the Java Runtime Environment, which may prevent conflicts with Java versions that may already be installed on your machine.<br>

* **Tracer:** Just like BEAST2, Tracer is written in Java and should work on your system without problems. The latest version of the program can be downloaded for Mac OS X, Linux, or Windows from [http://tree.bio.ed.ac.uk/software/tracer/](http://tree.bio.ed.ac.uk/software/tracer/). The file with the extension `.dmg` is for Mac OS X, the one with the extension `.tgz` is for Linux, and the Windows version is the file ending in `.zip`.

* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by Andrew Rambaut is a very intuitive and useful tool for the visualization and (to a limited extent) manipulation of phylogenies encoded in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format. Executables for Mac OS X, Linux, and Windows are provided on [http://tree.bio.ed.ac.uk/software/figtree/](http://tree.bio.ed.ac.uk/software/figtree/).


<a name="starbeast2"></a>
## Bayesian species-tree inference with StarBEAST2

In the part of the tutorial, we are going to use the multi-species-coalescent model implementation of StarBEAST2 ([Ogilvie et al. 2017](https://academic.oup.com/mbe/article/34/8/2101/3738283)) to estimating a time-calibrated species tree from the set of twelve alignments. If you're not familiar with the BEAST2 environment yet, you might want to check the tutorials on [Bayesian Phylogenetic Inference](../bayesian_phylogeny_inference/README.md) and [Phylogenetic Divergence-Time Estimation](../divergence_time_estimation/README.md) before going through this tutorial.

* Download the compressed directory [`09.tgz`](data/09.tgz) containing the 72 filtered alignments produced in tutorial [Ortholog Detection](../ortholog_detection/README.md).

* If your browser did not already uncompress the directory, do so with the following command:

		tar -xzf 09.tgz

* As analyses with StarBEAST2 are relatively computationally demanding, we are going to limit the dataset to those alignments that contain sequence information for all of the eleven cichlid species. This is not the case for most of the 72 alignments in directory `09` because sequences for many genes were apparently not included in the transcriptome assembly for *Ophthalmotilapia ventralis* generated by [Baldo et al. (2011)](https://academic.oup.com/gbe/article/doi/10.1093/gbe/evr047/583924). To remove these alignments without sequence information for *Ophthalmotilapia ventralis*, we can use the Ruby script [`filter_genes_by_missing_data.rb`](src/filter_genes_by_missing_data.rb). Like the Ruby scripts used in tutorial [Ortholog Detection](../ortholog_detection/README.md), the first two command-line arguments required by this script are the names of a directory with sequence alignments and the name of another directory to which the filtered alignments should be written. In addition, a third arguments is required for the number of sequences per alignment that may be completely missing. In our case, this number should be zero to ensure that all alignments contain at least partial sequence information for all species. Thus, run the script with the following command:

		ruby filter_genes_by_missing_data.rb 09 10 0

	**Question 1:** How many alignments are left after this filtering step? [(see answer)](#q1)

* Open BEAUti and click on "Manage Packages" in the "File" menu. This should open the BEAST2 Package Manager in a new window as shown in the next screenshot. Scroll to the bottom of the list of available packages, select "StarBEAST2", and click on "Install/Upgrade". Also make sure that the bModelTest package is installed, which should be the case if you did the tutorial on [Bayesian Phylogenetic Inference](../bayesian_phylogeny_inference/README.md) earlier (you can tell that a package is installed if a version number is listed for it in the second column of the table). The CladeAge package is not going to be needed for this tutorial, so you could uninstall it now if you installed it earlier.<p align="center"><img src="img/beauti1.png" alt="BEAUti" width="700"></p>

* You may notice that the interface of BEAUti apparently has not yet changed after installing the StarBEAST2 package. This is because we still need to load one of several BEAUti templates that are provided by StarBEAST2. If you hover with the mouse over "Template" in BEAUti's "File" menu, you'll see that several templates are now available that apparently are connected to StarBEAST2, as shown in the next screenshot.<p align="center"><img src="img/beauti2.png" alt="BEAUti" width="700"></p>The templates provided by StarBEAST2 are "SpeciesTreeUCED", "SpeciesTreeRLC", "StarBeast2", and "SpeciesTreeUCLN" (the template named "StarBeast" refers to the older version of StarBEAST that is installed by default). Of the four templates for StarBEAST2, the one named "StarBeast2" implements a strict-clock model, the two templates named "SpeciesTreeUCED" and "SpeciesTreeUCLN" implement relaxed-clock models with exponentially or lognormally distributed rate variation, respectively, and the template named "SpeciesTreeRLC" implements a random local clock. In contrast to the first version of StarBEAST (\*BEAST), rate variation is modelled in StarBEAST2 as a species trait that applies to all genes instead of being estimated individually for each gene; this is one of the reasons why StarBEAST2 is much faster than StarBEAST (\*BEAST). The different clock models and their implementations are described in detail in [Ogilvie et al. (2017)](https://academic.oup.com/mbe/article/34/8/2101/3738283).<br>Because we are here investigating relationships of rather closely related species with comparable lifestyle and habitat, we'll assume that differences in their evolutionary rates are negligible. This choice will also be the more convenient one for this tutorial as relaxed-clock models would be computationally more demanding. Thus, select the template named "StarBeast2" to use the strict-clock model. Note, however, that if the analysis would be for a publication, it would be worth also exploring models of rate variation.<br>After clicking on the "StarBeast2" template, you should see that the tabs "Taxon sets", "Gene Ploidy", and "Population Model" have been added to the BEAUti window, as shown in the next screenshot.
<p align="center"><img src="img/beauti3.png" alt="BEAUti" width="700"></p>

* Click on "Import Alignment" in the "File" menu and select the twelve alignments of our dataset. The window should then look as shown below.<p align="center"><img src="img/beauti4.png" alt="BEAUti" width="700"></p>

* Select all partitions and click on "Link Clock Models" in the menu bar above the partitions table. However, unlike in tutorial [Bayesian Phylogenetic Inference](../bayesian_phylogeny_inference/README.md), leave the tree models unlinked. In contrast to the analyses in this other tutorial, we will also not split each partition according to codon positions. For a full analysis it could be worthwile to do so, but to reduce computational requirements, we will here use only one partition per gene.

* Move on to the tab named "Taxon Sets". To allow the estimation of population sizes, StarBEAST2 analyses are usually performed with sequences of more than one individual per species, and the table in tab "Taxon Sets" then allows one to specify which individuals belong to which species. Here, however, our dataset includes only sequences from a single individual of each species (this means that with our dataset we are unable to estimate population sizes reliably but we will avoid this problem by using a fixed population size, as will be described below). Nevertheless, it is required that we specify a species name for each taxon in our phylogeny, and these species names must not be identical to the names of the individuals. A simple solution is to reuse the names assigned to individuals also for the species and add a prefix such as "_spc" to the end of each name, as shown in the screenshot below.<p align="center"><img src="img/beauti5.png" alt="BEAUti" width="700"></p>

* We'll ignore the "Tip Dates" tab as in tutorial [Bayesian Phylogenetic Inference](../bayesian_phylogeny_inference/README.md). Have a look at the "Gene Ploidy" tab instead. This is where you can specify the ploidy for each gene, which would have to be adjusted if we had mitochondrial markers. Given that all genes are from the nuclear genome (and assuming that none are from sex chromosomes), the default ploidy of 2 is correct; thus, don't change anything in this tab.<p align="center"><img src="img/beauti6.png" alt="BEAUti" width="700"></p>

* Move on to the tab named "Population Model". As you can see from the selection in the drop-down menu at the top of this window, the default setting for the population model is "Automatic Population Size Integration". This default is a good option when multiple individuals are used per species; however, because our dataset includes only a single individual for each species, we should fix the population size to a reasonable estimate instead. To do so, select "Constant Populations" from the drop-down menu. In the field to the right of "Population Sizes", just leave the default value of 1.0. Even though unintuitive, this value does not directly specify the effective population size. Instead, this value needs to be scaled by the number of generations per time unit. Given that we will use 1 million years as the time unit in our analysis (as in the other tutorials), and assuming a generation time of 3 years for cichlids ([Malinsky et al. 2015](http://science.sciencemag.org/content/350/6267/1493)), there are 333,333 generations per time unit. Thus the value of 1.0 specified for the population size in fact translates to an assumed effective population size of 333,333, which is comparable to the population sizes estimated for African cichlid fishes in [Meyer et al. (2017)](https://academic.oup.com/sysbio/article-abstract/66/4/531/2670093). Make sure to leave the checkbox for to the left of "estimate" at the right of the window unticked so that the specified population size is in fact fixed. The window should then look as in the screenshot below.
<p align="center"><img src="img/beauti7.png" alt="BEAUti" width="700"></p>

* Continue to the "Site Model" tab. As in the other tutorials, select the "BEAST Model Test" model from the first drop-down menu to average over a set of substitution models, and choose the "namedExtended" set of substitution models from the second drop-down menu. Make sure to set a tick in the checkbox next to "estimate" at the right of the window to estimate the mutation rate of the first partition relative to those of other partitions. The BEAUti window should then look as shown in the next screenshot.<p align="center"><img src="img/beauti8.png" alt="BEAUti" width="700"></p>

* Still in the "Site Model" tab, select all partitions from the list at the left of the window as shown below, and click "OK" to copy the settings for the first partition also to all other partitions.<p align="center"><img src="img/beauti9.png" alt="BEAUti" width="700"></p>

* Move on to the "Clock Model" tab. The "Strict Clock" model should already be selected in the drop-down menu; don't change this selection. However, we are going to estimate the rate of the clock through a calibration on the root of the phylogeny (this will be described below). To enable this estimation of the clock rate, set a tick next to "estimate" at the right of the window, as shown in the next screenshot.<p align="center"><img src="img/beauti10.png" alt="BEAUti" width="700"></p>

* Continue to the "Priors" tab. At the very top of the list, select the "Birth Death Model" from the first drop-down menu next to "Tree.t:Species" to allow the possibility that extinction occurred during the evolution of cichlid fishes, as shown in the screenshot below.<p align="center"><img src="img/beauti11.png" alt="BEAUti" width="700"></p>

* Then, scroll to the very bottom of the list in the "Priors" tab. To specify a time calibration, click the "+ Add Prior" button at the very bottom. This should open a pop-up window as shown in the next screenshot. Keep the selection of "MRCA prior" in the drop-down menu and click "OK".<p align="center"><img src="img/beauti12.png" alt="BEAUti" width="700"></p>

* This should open another pop-up window in which you should select a tree to which the calibration should be applied. Unfortunately, the window lists the names of all gene trees next to each other and thus may be wider than the computer screen. The button for the selection of the species tree, named "Tree.t:Species", is located at the very right of the window, as shown in the screenshot below. Thus, to select the species tree, you might have to move the window to the left until you can click on this button.<p align="center"><img src="img/beauti13.png" alt="BEAUti" width="700"></p>

* Once you clicked the button, yet another pop-up window should open in which you can specify the ingroup of the clade that is to be time calibrated. We are going to calibrate the divergence of Neotropical and African cichlid fishes according to the results obtained in tutorial [Phylogenetic Divergence-Time Estimation](../divergence_time_estimation/README.md). And since the dataset used in the current tutorial contains no species besides Neotropical and African cichlid fishes, their divergence represents the root of the phylogeny. Thus, to constrain the age of the root, select all species names from the left side of the pop-up window and click the `>>` button to move all of them into the ingroup on the right-hand side of the pop-up window. In the field at the top of the pop-up window next to "Taxon set label", specify "All" as the name of this group, as shown in the screenshot below. Then, click "OK".<p align="center"><img src="img/beauti14.png" alt="BEAUti" width="700"></p>

* In tutorial [Phylogenetic Divergence-Time Estimation](../divergence_time_estimation/README.md), the age of the divergence of Neotropical and African cichlid fishes was estimated as 65 Ma, with a confidence interval ranging from 55 to 75 Ma. To implement this age as a calibration for the current phylogeny, we can select a normally distributed prior density with a mean of 65 and a standard deviation of 5.1. To do so, choose "Normal" from the drop-down menu to the right of "All.prior" and write "65" and "5.1" in the fields to the right of "Mean" and "Sigma", as in the screenshot below.<p align="center"><img src="img/beauti15.png" alt="BEAUti" width="700"></p>

* Scroll to the right of the window to see the shape of the specified prior density for the age of the divergence of Neotropical and African cichlid fishes, as shown in the next screenshot.<p align="center"><img src="img/beauti16.png" alt="BEAUti" width="700"></p>As listed in the summary statistics below the density plot, the 2.5% quantile is at exactly 55.0 Ma and the 97% quantile is at 75.0 Ma; thus 95% of the prior mass lie between these values.

* Finally, continue to the "MCMC" tab. Depending on the time available for this analysis, set the chain length to 100 million as shown in the next screenshot or even to 1 billion. In the first case, the analysis will take about an hour but will not reach stationarity. Setting the chain length to 1 billion could make sense if you plan to run the analysis overnight. The chain will then approach stationarity but the analysis might still not be fully complete. Note that if you decide to run a shorter chain, you could use the results of my analysis (linked below) for the rest of the tutorial.<br>Set the value for "Store Every" to 50,000. Click on the black triangle to the left of "tracelog" to open the settings for the log file, and specify "starbeast.log" as the file name and 50,000 as the log frequency (to the right of "Log Every") as shown in the screenshot below.<p align="center"><img src="img/beauti17.png" alt="BEAUti" width="700"></p>

* Click on the next black triangle to see the settings for the species-tree log file. Set the name of this file to "starbeast_species.trees" and the log frequency again to 50,000 as in the next screenshot.<p align="center"><img src="img/beauti18.png" alt="BEAUti" width="700"></p>

* The settings for the log output to the screen ("screenlog") do not need to be changed. But click on the triangle below it to open the settings for the log file for the first gene tree, named "t:ENSDARG00000028507". Here, add "starbeast_" before the default file, and specify again a log frequency of 50,000 as shown in the screenshot below.<p align="center"><img src="img/beauti19.png" alt="BEAUti" width="700"></p>

* Repeat the above step also for the log files of all other gene trees.

* Then, save the file using "Save As" in BEAUti's "File" menu and name it "starbeast.xml".

* Before we analyze file `starbeast.xml` with BEAST2, we still need to make small adjustments to the file. This is because for some reason, BEAUti writes a parameter for the population size to the file and also places operators and a prior density on this parameter, even though we specified to use a constant population size. As a result, the prior probability would change throughout the MCMC analysis without any impact on the likelihood, which massively increases the number of iterations required to reach stationarity. Thus, we'll completely remove this parameter from the BEAST2 input file. To do so, open the file in a text editor find the following line on which the parameter is introduced, and delete the line:

		<parameter id="constPopMean.Species" lower="0.0" name="stateNode">1.0</parameter>
				
	The prior for the population-size parameter is specified with the three lines below. Search for these and also delete them.
	
		<prior id="constPopMeanPrior.Species" name="distribution" x="@constPopMean.Species">
			<OneOnX id="OneOnX.13" name="distr"/>
		</prior>
		
	The parameter is included in a multi-parameter operator named "updownAll:Species". Do not remove the whole operator, but only the following line in which the population-size parameter is included in it:
	
			<down idref="constPopMean.Species"/>
			
	A second operator acts exclusively on the population-size parameter. It is specified on the following line, which should be removed:
	
		<operator id="constPopMeanScale.Species" spec="ScaleOperator" parameter="@constPopMean.Species" scaleFactor="0.75" weight="1.0"/>
	
	Finally, the following line specifies that the population-size parameter should be included in the log output. Also remove this line:
	
		<log idref="constPopMean.Species"/>
		
	Save your changes to file `starbeast.xml`. This file is then ready to be analyzed with BEAST2.
		
* Open the program BEAST2, select file `starbeast.xml`, and click on "Run" to start the analysis.

As described above, this analysis may take between 1 and 10 hours depending on the length of the MCMC chain that you decided to specify. Instead of waiting for the analysis to finish, you could cancel it at some point and use the output of my analysis to complete this tutorial. However, you could in any case first continue with the next part of the tutorial which is independent of the results of the StarBEAST2 analysis, and you could keep the StarBEAST2 analysis running in the meantime.

[Bayesian species-tree inference with concatenation](#concatenation)
<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** After removing alignments that do not contain sequence information for all eleven cichlid species, twelve sequence alignments remain in the dataset, as can be shown with

		ls 10/*.fasta | wc -l