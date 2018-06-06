# Maximum-Likelihood Species-Tree Inference

A tutorial on maximum-likelihood species-tree inference based on gene trees

## Summary

Due to incomplete lineage sorting and recombination, different regions of the genomes of species may differ in their phylogenetic histories. Thus, a set of alignments of sequences from different genomic locations may support not just one but multiple different "gene trees" (this expression is commonly used regardless of whether the sequences actually represent genes or other types of markers). These potential differences among the true gene trees are usually ignored when phylogenetic inference is based on concatenation of multiple alignments. In contrast, differences among gene trees due to incomplete lineage sorting are explicitly accounted for in species-tree inference based on the multi-species-coalescent model. By accounting for these difference, programs implementing the multi-species coalescent model have been shown to be statistically consistent in the presence of incomplete lineage sorting, which is not the case for inference based on concatenation. Thus, phylogeny inference with the multi-species coalescent model is often more reliable than concatenation, particularly when rapidly speciating groups or taxa with large population sizes are investigated.

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Maximum-likelihood gene-tree inference with RAxML](#raxml)
* [Species-tree inference with ASTRAL](#astral)

<a name="outline"></a>
## Outline

In this tutorial, I will present how to quickly generate a set of gene trees based on maximum likelihood, using the software [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) ([Stamatakis 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)). In addition to the maximum-likelihood phylogeny, a set of bootstrap trees will also be generated for each gene. The sets of maximum-likelihood and bootstrap gene trees are then going to be used for species-tree inference under the multi-species-coalescent model, as implemented in the software [ASTRAL](https://github.com/smirarab/ASTRAL) ([Zhang et al. 2017](https://link.springer.com/chapter/10.1007%2F978-3-319-67979-2_4)).

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

* **RAxML:** [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) ([Stamatakis 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)) source code for Mac OS X and Linux, as well as precompiled executables for Windows, can be found on RAxML's github page [https://github.com/stamatak/standard-RAxML](https://github.com/stamatak/standard-RAxML). The installation of RAxML is also described in more detail in tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md).

* **ASTRAL:** The program [ASTRAL](https://github.com/smirarab/ASTRAL) ([Zhang et al. 2017](https://link.springer.com/chapter/10.1007%2F978-3-319-67979-2_4)) allows efficient and accurate estimation of the species tree based on a set of gene trees. The latest release of ASTRAL can be obtained from [https://github.com/smirarab/ASTRAL/releases](https://github.com/smirarab/ASTRAL/releases). Click on the link for "Source code (zip)" to download the full release including the source code as well as the compiled program as a Java jar file. Within the downloaded directory you'll find a zip file with the program name and its version number. Uncompress this file and open the uncompressed directory, there you should find the ASTRAL jar file, named `astral.5.5.6.jar` or similar. Rename this file so that it is simply named `astral.jar` and place it in a convenient location on your computer.

* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by Andrew Rambaut is a very intuitive and useful tool for the visualization and (to a limited extent) manipulation of phylogenies encoded in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format. Executables for Mac OS X, Linux, and Windows are provided on [http://tree.bio.ed.ac.uk/software/figtree/](http://tree.bio.ed.ac.uk/software/figtree/).

<a name="raxml"></a>
## Maximum-likelihood gene-tree inference with RAxML

As input for the species-tree analyses with ASTRAL, sets of gene trees are required. These gene trees could be generated with Bayesian approaches such as BEAST2; however, since ASTRAL anyway ignores branch lengths and uses gene-tree topologies only, time calibration of gene trees is not necessary. Thus, we may as well use RAxML to generate gene trees based on maximum likelihood, which is much faster than Bayesian approaches. If you have not used RAxML before, you may find a more detailed introduction to analyses with this program in tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md).

* Download the compressed directory [`09.tgz`](data/09.tgz) containing the 72 filtered alignments produced in tutorial [Ortholog Detection](../ortholog_detection/README.md).

* If your browser did not already uncompress the directory, do so with the following command:

		tar -xzf 09.tgz

* As noted above (and as you may have noticed if you ran tutorial [Ortholog Detection](../ortholog_detection/README.md)) almost all of these 72 gene alignments do not contain sequence information for *Ophthalmotilapia ventralis* ("ophven"), presumably because the transcriptome assembly generated by [Baldo et al. (2011)](https://academic.oup.com/gbe/article/doi/10.1093/gbe/evr047/583924) may have been incomplete. Instead of sequence information, these alignments therefore contain only missing data, coded with the gap symbol "-", for *Ophthalmotilapia ventralis*. Find out how many of the alignments do not contain any information for *Ophthalmotilapia ventralis* at all, using the following command:

		cat 09/*.fasta | grep -A 1 ophven | grep -v ophven | grep -e A -e C -e G -e T | wc -l
		
	**Question 1:** How many alignments contain information for *Ophthalmotilapia ventralis*? [(see answer)](#q1)

* Since RAxML can not infer the phylogeny from an alignment in which one or more sequences consist only of missing data, we will need to remove the "ophven" sequence from all those alignments in which it contains no information. This means that *Ophthalmotilapia ventralis* will only be included in a subset of the gene trees generated with RAxML; however, this will fortunately not be a problem for the species-tree inference because ASTRAL does not require that all gene trees contain the exact same set of taxa. To remove sequences that contain only missing information from all alignments, and at the same time translate all alignments into Phylip format, we can use the Python script [`convert.py`](src/convert.py) with the following command:

		for i in 09/*.fasta
		do
			gene_id=`basename ${i%.fasta}`
			python3 convert.py ${i} 09/${gene_id}.phy -f phylip -m 0.9
		done

	(the "-m" option specifies the maximally allowed proportion of missing data per sequence; setting it to 0.9 thus removes all sequences that are less than 10% complete).

* We can now use RAxML to generate maximum-likelihood gene trees for all alignments. To be able to later use bootstrapping with ASTRAL, we will also generate sets of trees from bootstrapped gene alignments in the same RAxML analysis, by using RAxML option "-f a" as in tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md). However, unlike in the other tutorial, we will now ensure that the same number of bootstrap replicates is conducted for each gene tree, and we set this number to 100 by specifying "-N 100". To start RAxML in a loop so that it analyzes one gene alignment after the other, use the following commands:

		for i in 09/*.phy
		do
			gene_id=`basename ${i%.phy}`
			raxml -s ${i} -n ${gene_id} -m GTRGAMMA -p ${RANDOM} -f a -x ${RANDOM} -N 100
		done

	(by setting the options "-p" and "-x" to "${RANDOM}", two random numbers will be generated and used as the seeds for the parsimony inference of the starting tree and for bootstrapping). The 72 RAxML analyses should finish within a few minutes.
	
* The RAxML analyses will have generated a large number of files containing run info, maximum-likelihood trees with and without node-support labels, as well as the set of trees from bootstrapped alignments of each gene. To clean up the directory and keep only the important files, use the following commands:

		rm RAxML_info.*
		rm RAxML_bipartitions.*
		rm RAxML_bipartitionsBranchLabels.*

* Next, have a look at one of the files containing the maximum-likelihood trees without node support, e.g. with the following command:

		less RAxML_bestTree.ENSDARG00000002952
	
	(then type `q` to return to the command line).<br>As you'll see, this file contains only a single line with the maximum-likelihood tree for gene ENSDARG00000002952 in Newick format.
	
* Since ASTRAL will require as input a single file containing all gene trees, combine all files with maximum-likelihood trees into a single file named `ml_best.trees`, using the following command:

		cat RAxML_bestTree.* > ml_best.trees

* To further clean up the directory, you could then also remove all files that contain the maximum-likelihood trees for single genes, using

		rm RAxML_bestTree.*

<a name="astral"></a>
## Species-Tree Inference with ASTRAL

ASTRAL infers the species tree by searching for the topology that agrees with the largest number of species quartets included in the set of gene trees. In addition, branch lengths are estimated in coalescent units based on the number of quartets that are in conflict with a branch of the species tree. Unlike concatenation, the species-tree approach of ASTRAL has been shown to be statistically consistent under the multi-species coalescent model ([Mirarab et al. 2014](https://academic.oup.com/bioinformatics/article/30/17/i541/200803)), meaning that its species-tree estimate is guaranteed to converge to the correct species tree with increasing size of the dataset. However, this statistical consistency is based on the assumption that gene-tree topologies are inferred without error, which may often not be the case. Thus, one might want to take uncertainty in the gene trees into account, which can be done with the sets of bootstrapped trees for each gene. In this case, ASTRAL uses the maximum-likelihood trees to infer the species-tree topology and branch lengths, and additional analyses of the bootstrapped trees are used to quantify node support as the proportions of these trees supporting a given node. However, newer versions of ASTRAL are also able to quantify node support completely without bootstrapped trees, based on the maximum-likelihood gene trees alone. Thus, uncertainty in the gene trees is then not taken into account at all. Nevertheless, according to the authors ([Sayyari and Mirarab 2016](https://academic.oup.com/mbe/article/33/7/1654/2579300)), simulations have shown that this approach works just as well as or even better than the first approach based on bootstrapped trees. Here, we will test both approaches for quantifying node support with ASTRAL.

* We'll first use the set of bootstrapped trees to estimate node support on the species tree. To do so, ASTRAL requires as input a single file with the names of all the files containing bootstrapped trees. We can generate such a file with the following command:

		ls RAxML_bootstrap.* > ml_boot.txt
		
* We can then run ASTRAL with two input files: The file containing the maximum-likelihood trees for each gene (`ml_best.trees`), and the file containing the names of all files with bootstrapped trees (`ml_boot.txt`). The first of these is to be specified with ASTRAL's option "-i", and the second should be given with option "-b". In addition, we'll use option "-o" to set the name of the output file to `species_boot.trees`:

		java -jar astral.jar -i ml_best.trees -b ml_boot.txt -o species_boot.trees
		
	ASTRAL should finish this analysis within a few seconds.

* Have a look at the output file [`species_boot.trees`](res/species_boot.trees) using a text editor (or again the command `less`). You'll see that it contains 102 lines. The first 100 of these lines represent species trees in Newick format estimated for each the 100 bootstrapped trees of each gene. On line 101 is a consensus tree for these 100 trees. Finally, the last line contains the species tree estimated from the maximum-likelihood gene trees, annotated with node support based on the first 100 trees. Thus, the last line contains the species tree that we'll use for interpretation.

* However, before visualizing the species tree in FigTree, first conduct the second ASTRAL analysis based on the maximum-likelihood trees alone. Do so using the following command:

		java -jar astral.jar -i ml_best.trees -o species_pp.tre
		
	The output file named [`species_pp.tre`](res/species_pp.tre) now contains just a single species tree, annotated with posterior probabilities as node support.
	
* For some reason (probably a bug), the file [`species_pp.tre`](res/species_pp.tre) may not open in FigTree. However, since we are interested only in the last of the trees from file [`species_pp.tre`](res/species_pp.tre) as well as the tree from file [`species_pp.tre`](res/species_pp.tre), we'll generate a new file named `species.trees` that contains both of these two trees using the following commands:
		
		tail -n 1 species_boot.trees > species.trees
		cat species_pp.tre >> species.trees
		
* Try opening file [`species.trees`](res/species.trees) in FigTree. If this should not work (the reason is a bug that does not allow decimal positions in node support in the way it is specified by ASTRAL), open the file in a text editor, copy its content, and paste it into an new empty FigTree window. Click on "Midpoint Root" in the "Tree" menu, orient the tree, and select "label" to display as node labels (see tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md) if you're unsure how to do this), which should then appear more or less as shown in the next screenshot. Note that depending on the random number seeds used for RAxML, different gene trees may have resulted in your analysis and may also have led to a species tree that differs from the one shown here.<p align="center"><img src="img/figtree1.png" alt="FigTree" width="600"></p>

* As you can see in the top left of the window next to "Current Tree:", the displayed tree is the first out of two. The one shown in the screenshot above is the tree with node-support values based on bootstrapping. To see the next tree, click on the "Next" button near the top right of the window. This should display the tree with node-support values based only on the maximum-likelihood gene trees, as shown in the next screenshot.<p align="center"><img src="img/figtree2.png" alt="FigTree" width="600"></p>Unsurprisingly, the phylogeny itself is identical between the two trees, because the species tree is in both cases inferred from the maximum-likelihood gene trees. Perhaps more surprisingly, the node-support values correlate very strongly between the two trees, indicating that ASTRAL in fact does not require bootstrapping to estimate node support. **Question 2:** Do the taxonomic groups of cichlids, as listed in the table at the beginning of this tutorial, appear monophyletic in the species tree estimated by ASTRAL? [(see answer)](#q2)



<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** As you'll see from the screen output of the command, only twelve alignments contain sequence information for *Ophthalmotilapia ventralis*.

<a name="q2"></a>

* **Question 2:** In fact, all tribes appear monophyletic in the ASTRAL species tree, and so do Neotropical cichlids (comprising *Amphilophus citrinellus*, "ampcit", and *Andinoacara coeruleopunctatus*, "andcoe") as well as African cichlids (all remaining species). However, the ASTRAL species tree does not answer whether or not *Ophthalmotilapia ventralis* ("ophven") is the sister group to Haplochromini plus Lamprologini, or whether it is closer to one of the two tribes. Adding further sequence data for *Ophthalmotilapia ventralis*, for which only twelve gene sequences were used for the analysis, would probably help to infer its position more reliably. The species tree further supports a nested position of Lake Malawi and Lake Victoria species (*Metriaclima zebra*, "metzeb", and *Pundamilia nyererei*, "punnye", respectively) within the Lake Tanganyika radiation, a pattern that has long been known ([Salzburger et al. 2005](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-5-17)).