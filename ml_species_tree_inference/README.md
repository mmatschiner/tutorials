# Maximum-Likelihood Species-Tree Inference

A tutorial on maximum-likelihood species-tree inference based on gene trees

## Summary

Due to incomplete lineage sorting and recombination, different regions of the genomes of species may differ in their phylogenetic histories. Thus, a set of alignments of sequences from different genomic locations may support not just one but multiple different "gene trees" (this expression is commonly used regardless of whether the sequences actually represent genes or other types of markers). These potential differences among the true gene trees are usually ignored when phylogenetic inference is based on concatenation of multiple alignments. In contrast, differences among gene trees due to incomplete lineage sorting are explicitly accounted for in species-tree inference based on the multi-species-coalescent model. By accounting for these difference, programs implementing the multi-species coalescent model have been shown to be statistically consistent in the presence of incomplete lineage sorting, which is not the case for inference based on concatenation. Thus, phylogeny inference with the multi-species coalescent model is often more reliable than concatenation, particularly when rapidly speciating groups or taxa with large population sizes are investigated.

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Maximum-likelihood gene-tree inference with IQ-TREE](#iqtree)
* [Species-tree inference with ASTRAL](#astral)
* [Concordance analyses with IQ-TREE](#concordance)

<a name="outline"></a>
## Outline

In this tutorial, I will present how to quickly generate a set of gene trees based on maximum likelihood, using the software [IQ-TREE](http://www.iqtree.org) ([Nguyen et al. 2015](https://academic.oup.com/mbe/article/32/1/268/2925592)). In addition to the maximum-likelihood phylogeny, a set of bootstrap trees will also be generated for each gene. The sets of maximum-likelihood and bootstrap gene trees will then be used jointly for species-tree inference under the multi-species-coalescent model, as implemented in the software [ASTRAL](https://github.com/smirarab/ASTRAL) ([Zhang et al. 2017](https://link.springer.com/chapter/10.1007%2F978-3-319-67979-2_4)). Finally, the proportion of genes and sites supporting the species tree will be quantified with concordance analyses, again using IQ-TREE.

<a name="dataset"></a>
## Dataset

The dataset used in this tutorial is the set of alignments for 72 genes produced in tutorial [Ortholog Detection](../ortholog_detection/README.md). In brief, this dataset includes sequences for eleven cichlid species, two of which represent Neotropical cichlids while the remaining nine species are from Africa. The focus of the taxon set is on cichlids of the rapid radiations in the East African lakes Tanganyika, Malawi, and Victoria. The table below lists all species included in the set of alignments. Note, however, that the sequence data for *Ophthalmotilapia ventralis* were extracted from a transcriptome assembly whereas genome assemblies were available for all other species.

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

* **IQ-TREE:** Precompiled binaries for Mac OS X, Linux, and Windows are available on [http://www.iqtree.org/#download](http://www.iqtree.org/#download). To install IQ-TREE on any of these systems, download the version for your operating system, and decompress this file on your machine if necessary. In the decompressed directory, you'll find a subdirectory named `bin` and inside of this subdirectory should be a file named `iqtree` or `iqtree.exe`. To easily access this executable from the command line, also place it in a directory that is included in your [PATH](https://en.wikipedia.org/wiki/PATH_(variable))).

* **ASTRAL:** The program [ASTRAL](https://github.com/smirarab/ASTRAL) ([Zhang et al. 2017](https://link.springer.com/chapter/10.1007%2F978-3-319-67979-2_4)) allows efficient and accurate estimation of the species tree based on a set of gene trees. The latest release of ASTRAL can be obtained from [https://github.com/smirarab/ASTRAL/releases](https://github.com/smirarab/ASTRAL/releases). Click on the link for "Source code (zip)" to download the full release including the source code as well as the compiled program as a Java jar file. Within the downloaded directory you'll find a zip file with the program name and its version number. Uncompress this file and open the uncompressed directory, there you should find the ASTRAL jar file, named `astral.5.5.6.jar` or similar. Rename this file so that it is simply named `astral.jar` and place it in a convenient location on your computer.

* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by Andrew Rambaut is a very intuitive and useful tool for the visualization and (to a limited extent) manipulation of phylogenies encoded in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format. Executables for Mac OS X, Linux, and Windows are provided on [https://github.com/rambaut/figtree/releases](https://github.com/rambaut/figtree/releases).

<a name="iqtree"></a>
## Maximum-likelihood gene-tree inference with IQ-TREE

As input for the species-tree analyses with ASTRAL, sets of gene trees are required. These gene trees could be generated with Bayesian approaches such as BEAST2; however, since ASTRAL anyway ignores branch lengths and uses gene-tree topologies only, time calibration of gene trees is not necessary. Thus, we may as well use IQ-TREE to generate gene trees based on maximum likelihood, which is much faster than Bayesian approaches. If you have not used IQ-TREE before, you may find a more detailed introduction to analyses with this program in tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md).

* Download the compressed directory [`09.tgz`](data/09.tgz) containing the 72 filtered alignments produced in tutorial [Ortholog Detection](../ortholog_detection/README.md).

* If your browser did not already uncompress the directory, do so with the following command:

		tar -xzf 09.tgz

* As noted above (and as you may have noticed if you ran tutorial [Ortholog Detection](../ortholog_detection/README.md)) almost all of these 72 gene alignments do not contain sequence information for *Ophthalmotilapia ventralis* ("ophven"), presumably because the transcriptome assembly generated by [Baldo et al. (2011)](https://academic.oup.com/gbe/article/doi/10.1093/gbe/evr047/583924) may have been incomplete. Instead of sequence information, these alignments therefore contain only missing data, coded with the gap symbol "-", for *Ophthalmotilapia ventralis*. Find out how many of the alignments do not contain any information for *Ophthalmotilapia ventralis* at all, using the following command:

		cat 09/*.fasta | grep -A 1 ophven | grep -v ophven | grep -e A -e C -e G -e T | wc -l
		
	**Question 1:** How many alignments contain information for *Ophthalmotilapia ventralis*? [(see answer)](#q1)

* Since IQ-TREE can not infer the phylogeny from an alignment in which one or more sequences consist only of missing data, we will need to remove the "ophven" sequence from all those alignments in which it contains no information. This means that *Ophthalmotilapia ventralis* will only be included in a subset of the gene trees generated with IQ-TREE; however, this will fortunately not be a problem for the species-tree inference because ASTRAL does not require that all gene trees contain the exact same set of taxa. To remove sequences that contain only missing information from all alignments, and at the same time translate all alignments into Nexus format, we can use the Python script [`convert.py`](src/convert.py) with the following command:

		for i in 09/*.fasta
		do
			gene_id=`basename ${i%.fasta}`
			python3 convert.py ${i} 09/${gene_id}.nex -f nexus -m 0.9
		done

	(the "-m" option specifies the maximally allowed proportion of missing data per sequence; setting it to 0.9 thus removes all sequences that are less than 10% complete).

* We can now use IQ-TREE to generate maximum-likelihood gene trees for all alignments. To be able to later use bootstrapping with ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option `-bb` as in tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md). Also as in this other tutorial, we do not specify a substition model with the `-m` option and therefore allow IQ-TREE to automatically select the best-fitting model. Unlike in the other tutorial, however, we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option `--wbt`. To start IQ-TREE in a loop so that it analyzes one gene alignment after the other, use the following commands:

		for i in 09/*.nex
		do
			iqtree -s ${i} -bb 1000 --wbt
		done

	The 72 IQ-TREE analyses should finish within a few minutes.
	
* The IQ-TREE analyses will have generated a number of files for each gene, containing run info, a distance matrix, starting trees, and so on. The only output files required for our further analyses are those ending in `.treefile` (the maximum-likelihood gene trees with branch lengths) and `.ufboot` (the set of bootstrap trees without branch lengths). To clean up the directory and keep only the important files, use the following commands:

		rm 09/*.bionj
		rm 09/*.ckp.gz
		rm 09/*.contree
		rm 09/*.iqtree
		rm 09/*.log
		rm 09/*.mldist
		rm 09/*.model.gz
		rm 09/*.splits.nex
		rm 09/*.uniqueseq.phy		

* Next, have a look at one of the files containing the maximum-likelihood trees, e.g. with the following command:

		less 09/ENSDARG00000002952.nex.treefile
	
	(then type `q` to return to the command line).<br>As you'll see, this file contains only a single line with the maximum-likelihood tree for gene ENSDARG00000002952 in Newick format.
	
* Since ASTRAL will require as input a single file containing all gene trees, combine all files with maximum-likelihood trees into a single file named `ml_best.trees`, using the following command:

		cat 09/*.treefile > ml_best.trees

* To further clean up the directory, you could then also remove all files that contain the maximum-likelihood trees for single genes, using

		rm 09/*.treefile

<a name="astral"></a>
## Species-Tree Inference with ASTRAL

ASTRAL infers the species tree by searching for the topology that agrees with the largest number of species quartets included in the set of gene trees. In addition, branch lengths are estimated in coalescent units based on the number of quartets that are in conflict with a branch of the species tree. Unlike concatenation, the species-tree approach of ASTRAL has been shown to be statistically consistent under the multi-species coalescent model ([Mirarab et al. 2014](https://academic.oup.com/bioinformatics/article/30/17/i541/200803)), meaning that its species-tree estimate is guaranteed to converge to the correct species tree with increasing size of the dataset. However, this statistical consistency is based on the assumption that gene-tree topologies are inferred without error, which may often not be the case. Thus, one might want to take uncertainty in the gene trees into account, which can be done with the sets of bootstrap trees for each gene. In this case, ASTRAL uses the maximum-likelihood trees to infer the species-tree topology and branch lengths, and additional analyses of the bootstrap trees are used to quantify node support as the proportions of these trees supporting a given node. However, newer versions of ASTRAL are also able to quantify node support completely without bootstrapped trees, based on the maximum-likelihood gene trees alone. Thus, uncertainty in the gene trees is then not taken into account at all. Nevertheless, according to the authors ([Sayyari and Mirarab 2016](https://academic.oup.com/mbe/article/33/7/1654/2579300)), simulations have shown that this approach works just as well as or even better than the first approach based on bootstrapped trees. Here, we will test both approaches for quantifying node support with ASTRAL.

* We'll first use the set of bootstrapped trees to estimate node support on the species tree. To do so, ASTRAL requires as input a single file with the names of all the files containing bootstrapped trees. We can generate such a file with the following command:

		ls 09/*.ufboot > ml_boot.txt
		
* We can then run ASTRAL with two input files: The file containing the maximum-likelihood trees for each gene (`ml_best.trees`), and the file containing the names of all files with bootstrapped trees (`ml_boot.txt`). The first of these is to be specified with ASTRAL's option "-i", and the second should be given with option "-b". In addition, we'll use option "-o" to set the name of the output file to `species_boot.trees`:

		java -jar astral.jar -i ml_best.trees -b ml_boot.txt -o species_boot.trees
		
	ASTRAL should finish this analysis within a few seconds.

* Have a look at the output file [`species_boot.trees`](res/species_boot.trees) using a text editor (or again the command `less`). You'll see that it contains 102 lines. The first 100 of these lines represent species trees in Newick format estimated for each the first 100 bootstrapped trees of each gene (by default, ASTRAL uses only these; if we wanted to use more than 100 bootstrap trees per gene we would need to use the `-r` option, e.g. `-r 1000`). On line 101 is a consensus tree for these 100 trees. Finally, the last line contains the species tree estimated from the maximum-likelihood gene trees, annotated with node support based on the bootstrap trees. Thus, the last line contains the species tree that we'll use for interpretation.

* However, before visualizing the species tree in FigTree, first conduct the second ASTRAL analysis based on the maximum-likelihood trees alone. Do so using the following command:

		java -jar astral.jar -i ml_best.trees -o species_pp.tre
		
	The output file named [`species_pp.tre`](res/species_pp.tre) now contains just a single species tree, annotated with posterior probabilities as node support.
	
* Since we are interested only in the last of the trees from file [`species_pp.tre`](res/species_pp.tre) as well as the tree from file [`species_pp.tre`](res/species_pp.tre), we'll generate a new file named `species.trees` that contains both of these two trees using the following commands:
		
		tail -n 1 species_boot.trees > species.trees
		cat species_pp.tre >> species.trees
		
* Open the file [`species.trees`](res/species.trees) in FigTree. Click on "Decreasing Node Order" in the "Tree" menu to orient the tree, and select "label" to display as node labels (see tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md) if you're unsure how to do this), which should then appear more or less as shown in the next screenshot. Note that depending on the randomly chosen starting points of the IQ-TREE analyses, slightly different gene trees may have resulted in your analysis and may also have led to a species tree that is not exactly identical to the one shown here.<p align="center"><img src="img/figtree1.png" alt="FigTree" width="600"></p>

* As you can see in the top left of the window next to "Current Tree:", the displayed tree is the first out of two. The one shown in the screenshot above is the tree with node-support values based on bootstrapping. To see the next tree, click on the "Next" button near the top right of the window. This should display the tree with node-support values based only on the maximum-likelihood gene trees, as shown in the next screenshot.<p align="center"><img src="img/figtree2.png" alt="FigTree" width="600"></p>Unsurprisingly, the phylogeny itself is identical between the two trees, because the species tree is in both cases inferred from the maximum-likelihood gene trees. Perhaps more surprisingly, the node-support values correlate very strongly between the two trees, indicating that ASTRAL in fact does not require bootstrapping to estimate node support. **Question 2:** Do the taxonomic groups of cichlids, as listed in the table at the beginning of this tutorial, appear monophyletic in the species tree estimated by ASTRAL? [(see answer)](#q2)


<a name="concordance"></a>
## Concordance analyses with IQ-TREE

The latest version of IQ-TREE (from v.1.7) ([Minh et al. 2018](https://www.biorxiv.org/content/10.1101/487801v1)) has a useful feature that allows the quick calculation of concordance factors, indicating how well each gene, or each site of the alignments, agrees with a certain topology. A very helpful description of these features can be found on [Robert Lanfear's blog](http://robertlanfear.com/blog/files/concordance_factors.html). In brief, gene-concordance factors (gCF) quantify the proportion of trees of a certain set of gene trees that contain a particular node of a user-provided reference tree. Similarly, site-concordance factors (sCF) indicate the proportion of sites in an alignment that support each node of the reference tree (specifically, at each site the support for a node is compared with the support for two alternative topologies, and the site is said to support the node if it agrees better with the current topology than with the two alternatives).

* Let's first check how many of the maximum-likelihood gene trees support each individual node in the species tree produced by ASTRAL. This can be done with the following command, in which we specify the reference tree with `-t`, provide a prefix for output files with `--prefix`, and activate the analysis of gene-concordance factors with the `--gcf` option, followed by the name of the file containing gene trees:

		iqtree -t species_pp.tre --gcf ml_best.trees --prefix gene_concordance

	The IQ-TREE output should report that a tree with concordance factors was written to file `gene_concordance.cf.tree`.
	
* Open file [`gene_concordance.cf.tree`](res/gene_concordance.cf.tree) in FigTree, orient the tree, and select "label" to display as node labels. The tree should then appear as shown in the next screenshot.<p align="center"><img src="img/figtree3.png" alt="FigTree" width="600"></p>

	You'll see that two types of support values are shown for each node, separated by forward slashes. The first of these values are familiar, they are the posterior probabilities that were calculated by ASTRAL when we ran it with the maximum-likelihood gene trees alone instead of using the bootstrap trees. Thus, these values do not result from the concordance analysis with IQ-TREE, but they are shown because IQ-TREE found them in the reference-tree file `species_pp.tre`. The values after the forward slashes, however, are the gene-concordance factors calculated by IQ-TREE from the set of gene trees in file `ml_best.trees`. As you can see the gCF values appear to correlate with the ASTRAL's posterior probabilities but are generally lower. For example, it appears that only about 40% of the gene trees support the sister-group relationship between *Neolamprologus brichardi* (neobri) and *Neolamprologus olivaceous* (neooli) even though ASTRAL estimated the posterior probability of this relationship close to 1. This is not necessarily a contradiction, though, because incomplete lineage sorting, which is accounted for by ASTRAL, is expected to lead to disagreement between the species tree and gene trees particularly for rapidly diverging species.
	
* To also calculate site-concordance factors, we'll first need to generate a concatenated version of all gene alignments. We can do that with the Ruby script [`concatenate.rb`](src/concatenate.rb), using the following command:

		ruby concatenate.rb -i 09/*.nex -o concatenated.nex -f nexus

	(`-i` specifies the set of gene alignments, `-o` specifies the name of the concatenated output file, and `-f` specifies the format of the output file).
	
* Have a look at the concatenated file [`concatenated.nex`](res/concatenated.nex) using the `less` command:

		less concatenated.nex
		
	(type `q` to return to the command line). **Question 3:** How many sites are now included in the concatenated alignment? [(see answer)](#q3)
	
* By specifying the concatenated alignment with the `-s` option and by activating the calculation of site-concordance factors with the `--scf` option, followed by a maximum number of species quartets to be considered in the analysis (100 is by far sufficient), we can now calculate gene and site-concordance factors in the same IQ-TREE analysis with the following command:

		iqtree -t species_pp.tre --gcf ml_best.trees -s concatenated.nex --scf 100 --prefix gene_and_site_concordance

* Open file [`gene_and_site_concordance.cf.tree`](res/gene_and_site_concordance.cf.tree) in FigTree, and again orient the tree and select "label" to display as node labels. The tree should then appear as shown in the next screenshot.<p align="center"><img src="img/figtree4.png" alt="FigTree" width="600"></p>

	You'll see that three different support values are now shown for each node, these are the posterior probabilities from ASTRAL, the gene-concordance factors, and the site-concordance factors.  **Question 4:** Which proportion of sites supports the monophyly of Lamprologini? [(see answer)](#q4)

<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** As you'll see from the screen output of the command, only twelve alignments contain sequence information for *Ophthalmotilapia ventralis*.

<a name="q2"></a>

* **Question 2:** In fact, all tribes appear monophyletic in the ASTRAL species tree, and so do Neotropical cichlids (comprising *Amphilophus citrinellus*, "ampcit", and *Andinoacara coeruleopunctatus*, "andcoe") as well as African cichlids (all remaining species). However, the ASTRAL species tree does not answer whether or not *Ophthalmotilapia ventralis* ("ophven") is the sister group to Haplochromini plus Lamprologini, or whether it is closer to one of the two tribes. Adding further sequence data for *Ophthalmotilapia ventralis*, for which only twelve gene sequences were used for the analysis, would probably help to infer its position more reliably. The species tree further supports a nested position of Lake Malawi and Lake Victoria species (*Metriaclima zebra*, "metzeb", and *Pundamilia nyererei*, "punnye", respectively) within the Lake Tanganyika radiation, a pattern that has long been known ([Salzburger et al. 2005](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-5-17)).

<a name="q3"></a>

* **Question 3:** The number following the "nchar" keyword on the fourth line indicates that 73623 sites are included in the concatenated alignment.

<a name="q4"></a>

* **Question 4:** As shown in the table at the beginning of this tutorial, the tribe called Lamprologini is in our dataset represented by four species: *Neolamprologus brichardi* (neobri), *Neolamprologus marunguensis* (neomar), *Neolamprologus gracilis* (neogra), and *Neolamprologus olivaceous* (neooli). The most ancestral node among the four species has the label "1/90.3/92.9", meaning that 90.3% of the gene trees (65 of the 72 trees) and 92.9% of the alignment sites support this node.