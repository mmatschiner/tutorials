# Maximum-Likelihood Phylogeny Inference

A tutorial on phylogeny inference with maximum likelihood.

## Summary

As the name indicates, maximum-likelihood phylogeny inference aims to find the parameters of an evolutionary model that maximize the likelihood of observing the dataset at hand. The model parameters include the tree topology and its branch lengths but also all parameter of the substitution model (such as HKY or GTR) assumed in the inference. As the search space for these parameters is enormous when the dataset contains more than just a handful of taxa, all modern programs for maximum-likelihood phylogenetic inference apply heuristics to reach the maximum-likelihood parameter combination.

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Maximum-likelihood phylogeny inference with RAxML](#raxml)
* [Reading and visualizing tree files](#figtree)
* [Assessing node support with bootstrapping](bootstrap)


<a name="outline"></a>
## Outline

In this tutorial, I will present maximum-likelihood inference with one of the fastest programs developed for this type of analysis, the program [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) ([Stamatakis 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)). I will demonstrate how the reliability of nodes in the phylogeny can be assessed with bootstrapping ([Felsenstein 1985](https://www.jstor.org/stable/2408678)), how different substitution models can be applied to separate partitions, and how alignments of multiple genes can be concatenated to be jointly used in the same phylogenetic analysis.

<a name="dataset"></a>
## Dataset

The data used in this tutorial are the filtered versions of the alignments generated for 16s and rag1 sequences in tutorial [Multiple Sequence Alignment](../multiple_sequence_alignment/README.md). More information on the origin of the dataset can be found in that tutorial. RAxML requires input files in Phylip format, thus we will use the files [`16s_filtered.phy`](data/16s_filtered.phy) and [`rag1_filtered.phy`](data/rag1_filtered.phy).

<a name="requirements"></a>
## Requirements

* **RAxML:** Source code for Mac OS X and Linux, as well as precompiled executables for Windows, can be found on RAxML's github page [https://github.com/stamatak/standard-RAxML](https://github.com/stamatak/standard-RAxML). To install RAxML on any of these systems, download the [latest release](https://github.com/stamatak/standard-RAxML/releases), either in its zip or tar.gz-compressed version. Decompress this file on your machine.<br>
For installation on Linux, instructions are provided in the `README` file that you will find in the decompressed RAxML package. To compile the parallelized PTHREADS version of RAxML, try running

		make -f Makefile.AVX.PTHREADS.gcc

	If this should not work because your computer does not support [AVX](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions), you could try
	
		make -f Makefile.SSE3.PTHREADS.gcc

	If this also fails because your computer also doesn't support [SSE](https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions), then you should nevertheless be able to compile the slower standard version with
	
		make -f Makefile.PTHREADS.gcc

	On Mac OS X, you'll have to use the versions of the `Makefile` that end in `.mac`. This is not described in RAxML's `README` file. So to compile the sequential version, you'll have to run
	
		make -f Makefile.AVX.PTHREADS.mac
		
	or
	
		make -f Makefile.SSE3.PTHREADS.mac
		
	(no Mac version without AVX or SSE3 seems to be included in the latest releases).
		
	On Windows, just use the newest of the precompiled executables that you will find in a directory named `WindowsExecutables_v8.2.10` or similar, within the decompressed RAxML package.<br>

	The commands given in this tutorial will assume that you compiled the parallelized PTHREADS version of RAxML (rather than the sequential or the MPI version), that you named the file simply `raxml`, and that you placed it somewhere on your computer where your system can find it (i.e. in a directory that is included in your [PATH](https://en.wikipedia.org/wiki/PATH_(variable))). One way to guarantee this on Mac OS X or Linux is to place the executable in `/usr/local/bin`, for example using (if you compiled the AVX version)
	
		mv raxmlHPC-PTHREADS-AVX /usr/local/bin/raxml
		
	To verify that the RAxML executable can be found by your system, type
	
		which raxml
		
	If this command outputs a path such as `/usr/local/bin/raxml`, the executable can be found. As another check if RAxML is working as it should, type
	
		raxml -v
		
	and you should see the version number as well as a list of contributing developers. If you do, you're ready to start the tutorial.
	
* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by Andrew Rambaut is a very intuitive and useful tool for the visualization and (to a limited extent) manipulation of phylogenies encoded in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format. Executables for Mac OS X, Linux, and Windows are provided on [http://tree.bio.ed.ac.uk/software/figtree/](http://tree.bio.ed.ac.uk/software/figtree/).
	
<a name="raxml"></a>
## Maximum-likelihood phylogeny inference with RAxML

We will first generate a simple maximum-likelihood phylogeny only for the filtered 16s sequence alignment.

* To get an impression of the many options available in RAxML, have a look at the impressively long help text of the program:

		raxml -h
		
* Scroll back up to the beginning of the RAxML help text. Close to the top, you'll see that raxml could be started as easily as

		raxmlHPC -s sequenceFileName -n outputFileName -m substitutionModel 
	where "sequenceFileName" and "outputFileName" would have to be replaced with the actual sequence and output file names, and a substitution model would have to be chosen to replace "substitutionModel". Note that in our case, we would also start the program with "raxml", not "raxmlHPC", only because we named it that way.

* So, let's try to run a maximum-likelihood search, first for the 16S sequence data, using the alignment file in Phylip format [`16s_filtered.phy`](data/16s_filtered.phy). We'll start with as little command-line options as possible, and learn along the way which other options we need. I suggest doing so only for the reason that I find this easier than remembering all necessary options before the run. We'll use the GTRGAMMA model (the GTR model with gamma-distributed rate variation, as suggested for the 16S alignment by the model selection done in tutorial [`substitution_model_selection`](../substitution_model_selection/README.md)), and we choose `16s_filtered.out` as part of all result file names:

		raxml -s 16s_filtered.phy -n 16s_filtered.out -m GTRGAMMA 
* As you'll see, this minimalistic choice of options does not seem to be sufficient, and RAxML asks us to specify a random number seed with the option "-p". Before doing so, make sure to remove the log file (`RAxML_info.16s_filtered.out`) that RAxML just wrote to the current directory:

		rm RAxML_info.16s_filtered.out

	Then, try running RAxML again, this time with a random number seed:

		raxml -s 16s_filtered.phy -n 16s_filtered.out -m GTRGAMMA -p 123
		
* RAxML should finish the analysis within a few seconds and present output as shown in the screenshot below. **Question 1:** Are the inferred stationary frequencies of the four nucleotides (here called the "base frequencies") comparable to those inferred with PAUP\* in tutorial [`substitution_model_selection`](../substitution_model_selection/README.md))? [(see answer)](#q1) As you can see from RAxML's output, the best-scoring maximum-likelihood tree was written to file `RAxML_bestTree.16s_filtered.out`.
<p align="center"><img src="img/raxml1.png" alt="RAxML" width="600"></p>

<a name="figtree"></a>
## Reading and visualizing tree files

In this part of the tutorial, we will explore how phylogenetic trees are encoded in Newick format, the format used by almost all phylogenetic sofware, and we will visualize the maximum-likelihood phylogeny generated with RAxML with the program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). Fun fact: The Newick format is named after the [Newick's restaurant](http://www.newicks.com) in Dover, New Hampshire, where Joe Felsenstein and other developers of the format ["enjoyed the meal of lobsters"](http://evolution.genetics.washington.edu/phylip/newicktree.html) in 1986. A good explanation of the format and information on its origin can be also be found [here](http://evolution.genetics.washington.edu/phylip/newicktree.html).

* Open the file [`RAxML_bestTree.16s_filtered.out`](res/RAxML_bestTree.16s_filtered.out) in a text editor, or on the command line using for example the `less` command:

		less RAxML_bestTree.16s_filtered.out

	You'll see a long string containing the taxon IDs, each of which is followed by a codon and a number, and together with these, the taxon IDs are embedded in parentheses. As an example, a short (and simplified) segment of the string is this:
	
		(Mugilxxcephalu:0.16,(Synbranmarmora:0.29,Ambassispcxxxx:0.05):0.01)

* Open the program FigTree, copy the above short part of the tree string, and paste it into the new FigTree window. You'll see a phylogeny of the three taxa *Mugil cephalus* ("Mugilxxcephalu"), *Synbranchus marmoratus* ("Synbranmarmora"), and *Ambassis* sp. ("Ambassispcxxxx"), as shown in the screenshot below.
<p align="center"><img src="img/figtree1.png" alt="FigTree" width="600"></p>

* For better visualization, increase the font size for tip labels in the panel at the left (click on the triangle to the left of "Tip Labels" to open it), and untick the checkbox next to "Scale Bar", as shown in the next screenshot.
<p align="center"><img src="img/figtree2.png" alt="FigTree" width="600"></p>

* Also tick the checkbox next to "Branch Labels" to display the branch lengths as in the next screenshot. Note that the label of the branch leading to *Mugil cephalus* ("Mugilxxcephalu") is not displayed due to a minor bug; it is hidden by the menu bar at the top. However, you'll recognize that the other branch lengths correspond to the numbers specified in the tree string after the colons:
 
		(Mugilxxcephalu:0.16,(Synbranmarmora:0.29,Ambassispcxxxx:0.05):0.01)
<p align="center"><img src="img/figtree3.png" alt="FigTree" width="600"></p>

	By comparing the tree string and the visualization, you'll see that the parentheses encode the relationships among taxa. For example, the pair of parentheses around `Synbranmarmora:0.29,Ambassispcxxxx:0.05` specifies that the taxon names listed inside of it form one monophyletic clade. Thus, *Synbranchus marmoratus* ("Synbranmarmora") and *Ambassis* sp. ("Ambassispcxxxx") are defined as sister taxa that are closer to each other than either of them is to *Mugil cephalus* ("Mugilxxcephalu").
	
* Next, open the complete phylogeny generated by RAxML (file [`RAxML_bestTree.16s_filtered.out`](res/RAxML_bestTree.16s_filtered.out)) in FigTree. It should look as shown in the below screenshot.
<p align="center"><img src="img/figtree4.png" alt="FigTree" width="600"></p>
The way the phylogeny is rooted in the above screenshot is arbitrary, because we so far did not specify how the rooting should be done either in the RAxML analysis or in the visualization with FigTree. So this phylogeny provides no evidence that *Acantharchus pomotis* ("Acanthpomotis") is the sister to all the other teleost fishes (only that is is first among these in the alphabet).

* To correct the rooting of the phylogeny, we can specify an outgroup. From the taxonomy of the species included in this dataset, we know that zebrafish (*Danio rerio*; "Danioxxrerioxx") is a member of the clade named "Otomorpha" whereas all other species belong to the clade named "Euteleosteomorpha" ([Betancur-R. et al. 2017](https://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-017-0958-3)). Thus, the correct root of the phylogeny must lie between zebrafish and the other taxa. To specify zebrafish as an outgroup, click on the branch leading to "Danioxxrerioxx", as shown in the next screenshot.
<p align="center"><img src="img/figtree5.png" alt="FigTree" width="600"></p>

* Then, with that branch being selected, click on the "Reroot" icon with the yellow arrow in the menu bar. The phylogeny should then look as shown in the next screenshot.
<p align="center"><img src="img/figtree6.png" alt="FigTree" width="600"></p>

* As a final change, we could sort the taxa according to node order. To do so, click "Decreasing node order" in FigTree's "Tree" menu. This should move "Danioxxrerioxx" to the top of the plot:
<p align="center"><img src="img/figtree7.png" alt="FigTree" width="600"></p>
It is almost surprising how well this phylogeny resolves the correct relationships among the 41 taxa (which are known rather well from more extensive studies based on large molecular datasets as well as morphology). **Question 2:** Do cichlids appear monophyletic in this phylogeny (to answer this, you may need to look up the [table in the Multiple Sequence Alignment](../multiple_sequence_alignment/README.md) tutorial)? [(see answer)](#q2) **Question 3:** And are Neotropical cichlids (*Cichla temensis*, *Geophagus brasiliensis*, *Herichthys cyanoguttatus*) monophyletic? [(see answer)](#q3)

<a name="bootstrap"></a>
## Assessing node support with bootstrapping

XXX





<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** In the model selection done with PAUP\*, the estimates for the stationary frequencies under the GTR model with gamma-distrbuted rate variation were 0.314 (A), 0.232 (C), 0.211 (G), and 0.242 (T). In contrast, these frequencies were estimated by RAxML as 0.275 (A), 0.252 (C), 0.244 (G), and 0.230 (T). Thus, these estimates differ actually quite a lot. The reason for these differences could be the relatively small alignment size leading to stochastic variation, or the fact that the phylogeny was assumed in the model selection analysis with PAUP\*, but was co-estimated with the model parameters by RAxML.

<a name="q2"></a>

* **Question 2:** No, cichlids are not monophyletic in this phylogeny even though many of them cluster together. But e.g. the placement of the killifish *Aplocheilus panchax* ("Aplochepanchax") as the sister species to the Malagasy cichlid *Ptychochromis grandidieri* ("Ptychocgrandid") shows that cichlids do not appear as monophyletic in this phylogeny (in contrast to results from much more extensive studies that clearly show that cichlids are in fact a monophyletic group; [Matschiner et al. 2017](https://academic.oup.com/sysbio/article/66/1/3/2418030); [Betancur-R. et al. 2017](https://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-017-0958-3)).

<a name="q3"></a>

* **Question 3:** Yes, the three representatives of Neotropical cichlids form one monophyletic clade in this phylogeny.