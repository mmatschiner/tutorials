# Genome-Wide Genealogy Estimation with SNP Data

A tutorial on the estimation of genealogies across the entire genome with SNP data

## Summary

In sets of closely related species, phylogenetic approaches based on sequence alignments can be problematic because individual short alignments may not contain sufficient information for phylogenetic inference, but longer alignments are likely to contain many recombination breakpoints that can bias phylogenetic inference when they are not accounted for. On the other hand, SNP-based phylogenetic inference such as that supported by the software SNAPP (the focus of tutorial [Divergence-Time Estimation with SNP Data](../divergence_time_estimation_with_snp_data/README.md)), does not account for introgression and may therefore also produce misleading phylogenetic estimates when introgression is present. Thus, neither alignment-based phylogenetic inference, nor SNP-based inference with a tool like SNAPP may provide unbiased phylogenies of rapidly diverging clades in which introgression occurs. However, a possible solution to this problem may come from an entirely different type of methods has recently seen a lot of development. These methods use genome-wide SNPs to estimate the so-called ancestral recombination graph (ARG), which represents a collection of genealogies that are all linked to segments on the genome and separated from each other by inferred recombination breakpoints. Even when recombination breakpoints are very frequent, the segments in between them may not too short for reliable inference with ARG methods. The reason for this is that by jointly estimating the genealogies of all segments, each genealogy can also be informed by its neighboring genealogies, from which it should differ only in the position of a single branch due to recombination.

Until recently, methods for the inference of ARGs were far too computationally demanding to be applied to more than just a handful of species. This has changed dramatically with the release of the programs [tsinfer](https://github.com/tskit-dev/tsinfer) by [Kelleher et al. (2018)](https://www.biorxiv.org/content/10.1101/458067v1.abstract) and [Relate](https://myersgroup.github.io/relate/index.html) by [Speidel et al. (2019)](https://www.biorxiv.org/content/10.1101/550558v1), which can be applied even to thousands of samples. Of the two programs, we will here use Relate because tsinfer so far does not estimate divergence times, but Relate does. Common to both programs is that they do not consider species as something different from populations, and thus genetic exchange is always allowed among the units into which individuals are grouped. But as the distinction between populations and closely related species is often arbitrary anyway, this does not mean that Relate should not be applied to species-level datasets. As the program is extremely new, this tutorial should be considered work in progress and may only serve to understand the first basic steps of analyses with Relate.

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Genotype phasing](#phasing)
* [Estimating genome-wide genealogies with Relate](#relate)
* [Estimating population sizes and splits](#popsizes)


<a name="outline"></a>
## Outline

In this tutorial, I am going to present the estimation of genome-wide genealogies with the program Relate. The set of estimated genealogies will then serve to quantify population sizes and coalescence rates within and among species. Patterns of temporal variation in these coalescence rates may then allow to estimate relationships among species as well as their divergence times.


<a name="dataset"></a>
## Dataset

The SNP data used in this tutorial is the phased and imputed dataset generated in tutorial [Analysis of Introgression with Chromosome-Length Alignments](../analysis_of_introgression_with_chromosome_length_alignments/README.md). More detailed information about the origin of this dataset is given in the [Genotype phasing](../analysis_of_introgression_with_chromosome_length_alignments/README.md#phasing) section of this other tutorial. In brief, the dataset includes SNP data for the 28 samples of 14 cichlid species listed in the table below, and this data has already filtered based on read quality and depth. Only SNPs mapping to chromosome 5 of the tilapia genome assembly ([Conte et al. 2017](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3723-5)) are included in the dataset. In contrast to the analyses in tutorial [Analysis of Introgression with Chromosome-Length Alignments](../analysis_of_introgression_with_chromosome_length_alignments/README.md), we will not use the version of the dataset in which imputed genotypes were masked when they originally had missing data, because Relate does not allow sites with missing data.


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
| LJC9      | neocan     | *Neolamprologus cancellatus*  | Lamprologini  |
| LJD1      | neocan     | *Neolamprologus cancellatus*  | Lamprologini  |
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

* **Relate:** The software [Relate](https://myersgroup.github.io/relate/index.html) estimates genome-wide sets of genealogies and is extremely fast in doing so. Downloads for Mac OS X and Linux are provided on [https://myersgroup.github.io/relate/index.html](https://myersgroup.github.io/relate/index.html); however, installation on Windows is not supported. After downloading and decompressing the file containing the software, a collection of compiled programs (named `Relate`, `RelateMutationRate` etc.) can be found in the `bin` subdirectory inside of the decompressed directory.

* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) should already be installed if you followed the tutorials [Bayesian Phylogenetic Inference](../bayesian_phylogeny_inference/README.md), [Phylogenetic Divergence-Time Estimation](../divergence_time_estimation/README.md) or other tutorials. If not, you can download it for Mac OS X, Linux, and Windows from [https://github.com/rambaut/figtree/releases](https://github.com/rambaut/figtree/releases).

<!--XXX Replace with actual librarires XXX
* **Libraries for R:** Two R libraries will be required in this tutorial; these are [ape](https://cran.r-project.org/web/packages/ape/index.html) ([Paradis et al. 2004](https://academic.oup.com/bioinformatics/article/20/2/289/204981)) and [coda](https://cran.r-project.org/web/packages/coda/index.html) ([Plummer et al. 2006](http://oro.open.ac.uk/22547/)). If you're already familiar with R, just install these packages in the usual way. If not, the easiest way to do so might be via the command line. Type `R` to open the R environment interactively. Then, run the following commands:

		install.packages("ape", repos="http://ftp.gwdg.de/pub/misc/cran/", dependencies=T)
		install.packages("coda", repos="http://ftp.gwdg.de/pub/misc/cran/", dependencies=T)
		
	To ensure that both packages have been successfully installed, type these commands:
	
		library(ape)
		library(coda)

	If both commands result in no error messages, the two packages are ready to be used. Then, quit the R environment with `quit(save="no")`.
-->


<a name="relate"></a>
## Estimating genome-wide genealogies

Relate expects the input in the haplotype format used by the phasing program [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) ([Delaneau et al. 2012](https://www.nature.com/articles/nmeth.1785)), but Relate also comes with scripts that allow the conversion from VCF format into SHAPEIT's haplotype format. As these scripts, however, expect chromosomes to be named only with integers, we first need to replace the chromosome name "NC_031969" in the VCF file.

* To replace the chromosome name "NC_031969" in file [`NC_031969.f5.sub1.phased.vcf.gz`](../data/NC_031969.f5.sub1.phased.vcf.gz) with the number "5" (because the chromosome is in fact chromosome 5 of tilapia) and write a new file named `NC_031969.f5.sub1.phased.renamed.vcf.gz`, use the `sed` command:

		gunzip -c NC_031969.f5.sub1.phased.vcf.gz | sed "s/NC_031969/5/g" | gzip > NC_031969.f5.sub1.phased.renamed.vcf.gz

* To convert file `NC_031969.f5.sub1.phased.renamed.vcf.gz` into SHAPEIT's haplotype format, we can use the RelateFileFormats program that should be installed in the same directory as Relate itself. Run it with the following command:

		RelateFileFormats --mode ConvertFromVcf --haps NC_031969.f5.sub1.phased.renamed.haps --sample NC_031969.f5.sub1.phased.renamed.sample -i NC_031969.f5.sub1.phased.renamed

* In addition to the input files generated above, Relate is going to require a genetic map that quantifies recombination rates across the genome. A good map for our dataset, however, is currently not available. As a workaround, we will prepare a genetic map based on the assumption that the recombination rate is commpletely homogeneous across chromosome 5, and that 1 centimorgan corresponds to 1 Mbp. While this assumption is certainly not realistic, it was also made when phasing was performed with BEAGLE in tutorial [Analysis of Introgression with Chromosome-Length Alignments](../analysis_of_introgression_with_chromosome_length_alignments/README.md). To determine the robustness of Relate analyses to the specified genetic map, different maps could be tested; however, for the purpose of this tutorial, a map with homogeneous rates should be sufficient. The format of the map file is explained on the [Relate website](https://myersgroup.github.io/relate/input_data.html). Write the map with homogeneous rates using the text below, and save it to a new file named `flat.map`:

		pos COMBINED_rate Genetic_Map
		0       1     0
		38360022        0       38.360022

* In the inference of genome-wide genealogies with Relate, we also need to make assumptions about the effective population size and the mutation rate, but at least for the effictive population size, the specified value will be updated and estimated during the Relate analysis. As for the phasing in tutorial [Analysis of Introgression with Chromosome-Length Alignments](../analysis_of_introgression_with_chromosome_length_alignments/README.md), we can again use the population size estimate of 100,000, obtained in tutorial [Divergence-Time Estimation with SNP Data](../divergence_time_estimation_with_snp_data/README.md). Specifying an appropriate mutation rate is more difficult, because how well this rate corresponds to the dataset depends on how strictly the dataset is filtered. The good thing is that at least after the analysis, we will be able to judge how appropriate the mutation-rate estimate was, based on how well the resulting divergence times correspond to those that we estimated in tutorial [Divergence-Time Estimation with SNP Data](../divergence_time_estimation_with_snp_data/README.md). Tweaking the mutation rate and rerunning the analysis could then iteratively improve the results of the Relate analysis. As a start, we could use the mutation rate determined by sequencing father-mother-offspring trios from three different cichlid species by [Malinsky et al. (2018)](https://www.nature.com/articles/s41559-018-0717-x); this rate is 3.5 &times; 10<sup>-9</sup> per bp per generation. Keep in mind, however, that this rate is most likely an overestimate for our dataset because many true SNPs may have been removed in the filtering for quality and depth. Start the Relate analysis with the following command, in which the mutation rate of 3.5 &times; 10<sup>-9</sup> is specified with the `-m` option, the effective population size of 100,000 is specified with `-N`, and the genetic map file [`flat.map`](res/flat.map) is specified with `--map`:

		Relate --mode All -m 3.5e-9 -N 100000 --haps NC_031969.f5.sub1.phased.renamed.haps --sample NC_031969.f5.sub1.phased.renamed.sample --map flat.map -o NC_031969.f5.sub1.phased.renamed

	This analysis should finish within a few minutes.
	
* Check which files present in your current directory, for example using the `ls` command. You should see at least the following files:

		NC_031969.f5.sub1.phased.renamed.anc
		NC_031969.f5.sub1.phased.renamed.haps
		NC_031969.f5.sub1.phased.renamed.mut
		NC_031969.f5.sub1.phased.renamed.sample
		...
		
	The two files `NC_031969.f5.sub1.phased.renamed.anc` and `NC_031969.f5.sub1.phased.renamed.mut` are the output files of Relate. Of these, `NC_031969.f5.sub1.phased.renamed.anc` contains the inferred genealogies, and `NC_031969.f5.sub1.phased.renamed.mut` contains information about the mutations inferred to have occurred on those trees.
	
* Have a look at the content of `NC_031969.f5.sub1.phased.renamed.anc`, using, for example, the `less` command with option `-S` to avoid line wrapping. This should show more or less the following output:

		NUM_HAPLOTYPES 56
		NUM_TREES 38962
		0: 86:(1803.64969 2.077 0 65) 75:(4163.23310 6.758 0 575) 75:(4163.23310 13.454 0 575) 110:(
		65: 92:(2154.53825 0.000 65 86) 56:(4319.41568 6.758 0 575) 56:(4319.41568 13.454 0 575) 104
		86: 80:(3624.98532 1.500 86 272) 56:(5595.21322 6.758 0 575) 56:(5595.21322 13.454 0 575) 80
		111: 77:(3494.82820 1.500 86 272) 56:(4648.34200 6.758 0 575) 56:(4648.34200 13.454 0 575) 7
		120: 77:(3920.86808 1.500 86 272) 56:(4403.36588 6.758 0 575) 56:(4403.36588 13.454 0 575) 7
		137: 77:(3934.88503 1.500 86 272) 56:(4495.74078 6.758 0 575) 56:(4495.74078 13.454 0 575) 7
		177: 97:(3827.91020 1.500 86 272) 56:(4646.73315 6.758 0 575) 56:(4646.73315 13.454 0 575) 9
		188: 92:(4818.90426 1.500 86 272) 56:(5146.76792 6.758 0 575) 56:(5146.76792 13.454 0 575) 9
		196: 108:(4188.03328 1.500 86 272) 56:(4414.11292 6.758 0 575) 56:(4414.11292 13.454 0 575) 
		198: 108:(4474.93828 1.500 86 272) 56:(5706.95189 6.758 0 575) 56:(5706.95189 13.454 0 575) 
		201: 108:(3560.80238 1.500 86 272) 56:(5445.60224 6.758 0 575) 56:(5445.60224 13.454 0 575) 
		203: 108:(3698.91906 1.500 86 272) 56:(5081.65949 6.758 0 575) 56:(5081.65949 13.454 0 575) 
		211: 108:(4450.84729 1.500 86 272) 56:(4548.29507 6.758 0 575) 56:(4548.29507 13.454 0 575) 
		221: 107:(3818.56511 1.500 86 272) 56:(5379.28002 6.758 0 575) 56:(5379.28002 13.454 0 575) 
		240: 108:(3121.81418 1.500 86 272) 56:(4444.24124 6.758 0 575) 56:(4444.24124 13.454 0 575) 
		262: 97:(4155.74586 1.500 86 272) 56:(5258.00178 6.758 0 575) 56:(5258.00178 13.454 0 575) 9
		...
		
	As specified on the first two lines, this file includes 38,962 trees for 56 haplotypes. **Question 1:** Why are there 56 haplotypes? [(see answer)](#q1) Obviously, the format in which trees are saved in this file is not the Newick format.
	
* To visualize one randomly selected genealogy, we can use the RelateExtract program that should have been installed in the same directory as Relate itself. The following command extracts the tree that was inferred for chromosome position 1 Mbp and writes it to a new file in Newick format named `NC_031969.f5.sub1.phased.renamed_at_1000000.newick`:

		RelateExtract --mode TreeAtSNPAsNewick --anc NC_031969.f5.sub1.phased.renamed.anc --mut NC_031969.f5.sub1.phased.renamed.mut --bp_of_interest 1000000 -o NC_031969.f5.sub1.phased.renamed

* Open the tree file [`NC_031969.f5.sub1.phased.renamed_at_1000000.newick`](res/NC_031969.f5.sub1.phased.renamed_at_1000000.newick) in FigTree. The genealogy should look as shown in the next screenshot.<p align="center"><img src="img/figtree1.png" alt="FigTree" width="600"></p>

	As you can see, the tips in this genealogy are labeled with numbers instead of their sample IDs. It is possible, however, to link these numbers to sample IDs, because they correspond to the order in which samples are listed in file [`NC_031969.f5.sub1.phased.renamed.sample`](res/NC_031969.f5.sub1.phased.renamed.sample), and two consecutive numbers (e.g. "0" and "1") represent the two haplotypes of the same sample (in this case IZA1). We are going to translate the IDs later in this tutorial for a larger sample of trees. Another thing to notice in the genelogy shown in FigTree is the scale bar, which is in units of generations. The root age of this genealogy is around 259,000 generations, which is less than 1 million years when assuming a generation time of 3 years for cichlids as in tutorial [Divergence-Time Estimation with SNP Data](../divergence_time_estimation_with_snp_data/README.md). This is far younger than the divergence time estimated in this other tutorial, indicating that perhaps the mutation rate used for the Relate analysis was in fact too high for the dataset (there is no need to change it now, though).

* To figure out how the contents of files `NC_031969.f5.sub1.phased.renamed.anc` and `NC_031969.f5.sub1.phased.renamed.mut` correspond to each other, let's count the numbers of lines in both files:

		cat NC_031969.f5.sub1.phased.renamed.anc | wc -l
		cat NC_031969.f5.sub1.phased.renamed.mut | wc -l
		
	You'll see that file `NC_031969.f5.sub1.phased.renamed.mut` has many more lines than `NC_031969.f5.sub1.phased.renamed.anc`, 417,379 instead of 38,964.
	
* Now have a look at the content of file `NC_031969.f5.sub1.phased.renamed.mut`:

		less -S NC_031969.f5.sub1.phased.renamed.mut
		
	This should show something like the following content:
	
		snp;pos_of_snp;dist;rs-id;tree_index;branch_indices;is_not_mapping;is_flipped
		0;45158;13671;.;0;35 39;1;0;0;0;A/C;
		1;58829;29535;.;0;21 52 23 9;1;1;0;0;T/C;
		2;88364;7;.;0;;0;0;0;0;T/G;
		3;88371;9996;.;0;;0;0;0;0;C/T;
		4;98367;37;.;0;;0;0;0;0;T/A;
		5;98404;26364;.;0;26 38 24 36 50 49;1;0;0;0;G/C;
		6;124768;17088;.;0;9;0;0;0;1614.95;A/C;
		7;141856;40;.;0;2;0;1;0;4163.23;G/A;
		8;141896;12;.;0;6 4;1;0;0;0;G/T;
		9;141908;1760;.;0;;0;0;0;0;G/A;
		10;143668;6238;.;0;;0;0;0;0;A/C;
		11;149906;13;.;0;2 43 15;1;0;0;0;G/T;
		12;149919;6;.;0;26 24;1;0;0;0;C/T;
		13;149925;13898;.;0;3;0;0;0;8358.43;C/T;
		14;163823;245;.;0;;0;0;0;0;T/C;
		15;164068;16;.;0;;0;0;0;0;C/T;
		16;164084;27;.;0;109;0;0;5250.15;8358.43;G/A;
		17;164111;38645;.;0;;0;0;0;0;A/C;
		18;202756;87;.;0;44 39 10 69 31 12 14 9;1;0;0;0;C/T;
		...
		
	The rows in this file correspond to information from individual SNPs. It has twelve columns separated by semicolons, of which the second column specifies the position of the SNP on the chromosome. Importantly, the information in the fifth column specifies how SNPs are linked to the trees in file `NC_031969.f5.sub1.phased.renamed.anc`. All rows shown above contain "0" in the fifth column, indicating that they all represent SNPs associated with the first tree. You'll need to scroll down a bit (to line 66 in my case) to find the first SNP that is associated with the second tree, indicated by number "1" in the fifth column.
	

<a name="popsizes"></a>
## Estimating population sizes and splits

Once the genome-wide genealogies are inferred and mutations are mapped to them with Relate, the progam can use this wealth of information to estimate the sizes of different populations in the dataset. The estimation of split times, however, is only possible indirectly, through temporal variation in the coalescence rate between two populations. As mentioned at the beginning of this tutorial, Relate does not discriminate between species and populations and therefore there is not a single time point at which two groups of samples separate from each other, but instead the inferred temporal variation in coalescence rates more commonly indicates a gradual decline of among-population exchange, and the decision when exactly to call the two groups separate species may be arbitrary. Of course, this is in no way different with actual "species" in nature.

* To estimate population sizes and splits, we need to provide information about the population/species assignment to Relate. This is done with a "poplabels" file, the format of which is specified on the [Relate website](https://myersgroup.github.io/relate/input_data.html#.poplabels). Somewhat confusingly, two types of grouping are possible and expected, to indicate the "population" and the "group". In the example given on the Relate website, "GBR" (= Great Britain) is considered a "population" and "EUR" (= Europe) is a group. As these distinctions are somewhat arbitrary, we will only use the nominal species of the samples in our dataset for both the "population" and the "group". Furthermore, the sex of samples can be specified, and because we have that information for most samples, we can provide it. Males are coded with "1", females with "2", and samples without sex information can be marked with "NA". Thus, write the following information to a new file named `samples.poplabels`:

		sample	population	group	sex
		IZA1    astbur	astbur	1
		IZC5    astbur	astbur	2
		AUE7    altfas	altfas	1
		AXD5    altfas	altfas	2
		JBD5    telvit	telvit	1
		JBD6    telvit	telvit	2
		JUH9    neobri	neobri	2
		JUI1    neobri	neobri	1
		LJC9    neocan	neocan	NA
		LJD1    neocan	neocan	NA
		KHA7    neochi	neochi	1
		KHA9    neochi	neochi	2
		IVE8    neocra	neocra	1
		IVF1    neocra	neocra	2
		JWH1    neogra	neogra	2
		JWH2    neogra	neogra	1
		JWG8    neohel	neohel	2
		JWG9    neohel	neohel	1
		JWH3    neomar	neomar	1
		JWH4    neomar	neomar	2
		JWH5    neooli	neooli	2
		JWH6    neooli	neooli	1
		ISA6    neopul	neopul	1
		ISB3    neopul	neopul	2
		ISA8    neosav	neosav	1
		IYA4    neosav	neosav	2
		KFD2    neowal	neowal	2
		KFD4    neowal	neowal	1

* With the prepared "poplabels" file, we can then run the script `EstimatePopulationSize.sh` that should be located in the `scripts` directory of the Relate package. If for example, Relate itself is in `/Applications/Relate/bin/Relate` on your machine, then the script `EstimatePopulationSize.sh` should be placed in `/Applications/Relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh`. So, locate this script on your machine, and then execute the following command, in which we specify the prefix of the input files `NC_031969.f5.sub1.phased.renamed.anc` and `NC_031969.f5.sub1.phased.renamed.mut` with `-i`, the mutation rate of 3.5 &times; 10<sup>-9</sup> once again with `-m`, and the prefix of the output with `-o`. In addition, we need to specify the "poplabels" file [`samples.poplabels`](res/samples.poplabels) with `--poplabels` and the generation time of 3 years with `--years_per_gen`, otherwise the default of 28 years would be used. And finally, we'll use `--threshold` to specify a minimum number of 100 SNPs that should be associated with a genealogy, so that genealogies with less SNPs, which could be unreliable, are ignored in the analysis. As the default minimum SNP number per genealogy is the number of haplotypes, which is 56 in our case, increasing this threashold to 100 is going to reduce the number of genealogies used in the analysis and thus the run time of the script. The complete command to execute is then the following:

		PATH_TO_SCRIPT/EstimatePopulationSize.sh -i NC_031969.f5.sub1.phased.renamed -m 3.5e-9 --poplabels samples.poplabels --years_per_gen 3 --threshold 100 -o NC_031969.f5.sub1.phased.renamed.reestimated

	This script should take around 10 minutes to finish.
	
* The above command should have written a file named [`NC_031969.f5.sub1.phased.renamed.reestimated.coal`](res/NC_031969.f5.sub1.phased.renamed.reestimated.coal). Have a look at this file, for example with the `less` command and it `-S` option:

		less -S NC_031969.f5.sub1.phased.renamed.reestimated.coal
	
	You should something like the following content:
	
		altfas astbur neobri neocan neochi neocra neogra neohel neomar neooli neopul neosav 
		0 333.333 463.165 643.566 894.232 1242.53 1726.49 2398.95 3333.33 4631.65 6435.66 89
		0 0 0 0 2.57399e-06 4.62849e-06 1.47129e-06 5.73394e-07 6.22856e-07 2.50653e-05 4.15
		0 1 0 0 0 0 0 0 0 0 4.92089e-09 7.08465e-09 2.5501e-09 3.12072e-08 6.24374e-08 7.039
		0 2 0 0 0 0 0 0 0 0 3.44468e-08 1.41698e-08 3.18787e-08 1.30831e-07 1.66938e-07 2.72
		0 3 0 0 3.09884e-07 5.86596e-07 4.49828e-07 1.65144e-07 1.56979e-07 2.44754e-07 6.49
		0 4 0 0 0 0 0 0 0 0 0 2.83387e-08 3.5704e-08 1.07411e-07 2.20839e-07 2.83034e-07 4.2
		0 5 0 0 0 0 0 0 0 2.05218e-08 1.96839e-08 1.50554e-08 4.33556e-08 1.30382e-07 1.8216
		0 6 0 0 0 0 0 0 0 0 3.93676e-08 2.47971e-08 2.2953e-08 1.10169e-07 1.60652e-07 2.646
		0 7 0 0 0 0 0 0 0 3.42029e-08 2.46053e-08 1.417e-08 4.71819e-08 1.20277e-07 1.63636e
		0 8 0 0 0 0 0 0 0 0 1.96835e-08 2.12542e-08 5.61064e-08 1.78136e-07 1.79224e-07 2.85
		0 9 0 0 0 0 0 0 0 0 3.93679e-08 2.47972e-08 2.04028e-08 1.03738e-07 1.59654e-07 2.66
		0 10 0 0 0 0 0 0 0 4.10434e-08 1.96843e-08 1.59414e-08 4.0807e-08 1.15687e-07 1.6726
		0 11 0 0 0 0 0 0 0 0 4.92094e-08 2.1255e-08 2.86915e-08 1.48281e-07 1.93087e-07 3.04
		0 12 0 0 0 0 0 0 0 0 0 2.83387e-08 3.18786e-08 1.12002e-07 2.1787e-07 2.79706e-07 4.
		0 13 0 0 0 0 0 0 0 0 1.96835e-08 1.7712e-08 9.88328e-08 2.35569e-07 3.97802e-07 4.71
		1 0 0 0 0 0 0 0 0 0 4.92089e-09 7.08465e-09 2.5501e-09 3.12072e-08 6.24374e-08 7.039
		1 1 0 0 0 0 0 8.80312e-09 2.78919e-07 1.51536e-05 9.85812e-06 1.73112e-06 5.52349e-0
		1 2 0 0 0 0 0 0 0 6.84063e-09 0 7.08465e-09 2.5501e-09 1.65213e-08 7.59808e-08 6.991
		...
		
	In this file, the first two lines specify the names of the analyzed populations and the discrete time slices in numbers of generations within which coalescence rates have been estimated. In our analysis, the script `EstimatePopulationSize.sh` apparently used one time slice between the present and 333.3 generations ago, another time slice for the period between 333.3 generations ago and 463.2 generations ago, and so on. The breakpoints between these time slices are chosen automatically by the script and the number of time slices is by default set to 30, but we could have changed that number with the `--num_bins` option for `EstimatePopulationSize.sh`.

	The lines after the second line report in the first two columns two numbers corresponding to a population pair, and in the columns following these the estimated coalescence rates per time slice. The numbers used to indicate populations are according to the order in which population names were specified on the very first row of this file. So, for example, the fourth line of the [`NC_031969.f5.sub1.phased.renamed.reestimated.coal`](res/NC_031969.f5.sub1.phased.renamed.reestimated.coal) reports the coalescence rates for a pair of populations labeled on this line with "0" and "1". This means that the first and the second population from the list on the first line, "altfas" and "altbur" are compared. The coalescence rate is 0 in this pair in the first eight time slices, and 4.92089 &times; 10<sup>-9</sup> in the ninth time slice, and thus between 3333.3 generations ago and 4631.7 generations ago. The format of output files with the ending `.coal` is explained on the [Relate website](https://myersgroup.github.io/relate/modules.html#CoalescenceRate).
	
	The script `EstimatePopulationSize.sh` should also have written a file in PDF format named [`NC_031969.f5.sub1.phased.renamed.reestimated.pdf`](res/NC_031969.f5.sub1.phased.renamed.reestimated.pdf), in which changes in population sizes are plotted. However, the format of the plot, shown below, does not seem to be optimized for the numbers of populations in our dataset:<p align="center"><img src="img/NC_031969.f5.sub1.phased.renamed.reestimated.png" alt="FigTree" width="600"></p> The populations sizes shown in the above plot are calcualted as half of the inverse of the coalescence rate in a given time slice.
	
* As file [`NC_031969.f5.sub1.phased.renamed.reestimated.coal`](res/NC_031969.f5.sub1.phased.renamed.reestimated.coal) reports not only among-population coalescence rates but also within-population coalescence rates (for example in the third line, where population "0" is compared with population "0", itself), these plots can indeed inform about populations-size changes over time. For among-population comparisons, howver, it might be more meaningful to plot directly the inferred changes in coalescence rates, as these could tell us something about the divergence times between species. To prepare tables that can then be plotted in the R environment, use the Ruby script [`extract_coalescence_rates.rb`](src/extract_coalescence_rates.rb). This script expects the following four command-line arguments:

	* the name of the input file with ending `coal`.
	* the generation time,
	* the name of a species (population) for which coalescence times with all other species (populations) should be plotted,
	* the name of the output file to which a table will be written in a format that is more easily read with R.

	As a first test, we'll plot coalescence rates between *Astatotilapia burtoni* and all other species. As *Astatotilapia burtoni* is the outgroup species in our dataset, we expect that it has similar coalescence rates with all other species. Thus, execute the script with the following command:
	
		ruby extract_coalescence_rates.rb NC_031969.f5.sub1.phased.renamed.reestimated.coal 3 astbur coal_rates_astbur.txt
		
	(note that this script also replaces the breakpoints of time slices with the center of each time slice).
		
* Then, open the R environment, either using R Studio or another graphical user interface, or by typing `R` on the command line.

* Next, type the following commands in the R environment to read the table that we just wrote to file `coal_rates_astbur.txt` and produce a plot:

		table <- read.table("coal_rates_astbur.txt")
		pdf("coal_rates_astbur.pdf", height=7, width=7)
		plot(table$V1, table$V2, type="l", xlim=c(0,10000000), ylim=c(0,1E-5), xlab="Time", ylab="Coalescence rate", main="astbur")
		lines(table$V1, table$V3)
		lines(table$V1, table$V4)
		lines(table$V1, table$V5)
		lines(table$V1, table$V6)
		lines(table$V1, table$V7)
		lines(table$V1, table$V8)
		lines(table$V1, table$V9)
		lines(table$V1, table$V10)
		lines(table$V1, table$V11)
		lines(table$V1, table$V12)
		lines(table$V1, table$V13)
		lines(table$V1, table$V14)
		dev.off()
		
* Quit the R environment with `quit(save="no")`. The R commands should have written a plot named [`coal_rates_astbur.pdf`](res/coal_rates_astbur.pdf). This plot should look as shown below:<p align="center"><img src="img/coal_rates_astbur.png" alt="Relate" width="600"></p> As we can see, the coalescence rates between *Astatotilapia burtoni* and all other species in fact appear almost identical. There seems to be more variation with older ages, though, which could result from greater uncertainty in those estimates.

* For comparison, repeat the above steps with a species that is nested within the genus *Neolamprologus*, such as *N. olivaceous*. We would expect that the coalescence rates between *N. olivaceous* and other species of the dataset show much greater variation. Run the Ruby script `ruby extract_coalescence_rates.rb` again, but this time with `neooli` as the third argument, and with `coal_rates_neooli.txt` as the fourth argument:

		ruby extract_coalescence_rates.rb NC_031969.f5.sub1.phased.renamed.reestimated.coal 3 neooli coal_rates_neooli.txt

* Then, use the R environment again to produce a plot for *N. olivaceous*, using the following commands (unless you use R Studio or a similar tool, enter the R environment again by typing `R` and exit afterwards with `quit(save="no")`):

		table <- read.table("coal_rates_neooli.txt")
		pdf("coal_rates_neooli.pdf", height=7, width=7)
		plot(table$V1, table$V2, type="l", xlim=c(0,10000000), ylim=c(0,3E-5), xlab="Time", ylab="Coalescence rate", main="neooli")
		lines(table$V1, table$V3)
		lines(table$V1, table$V4)
		lines(table$V1, table$V5)
		lines(table$V1, table$V6)
		lines(table$V1, table$V7)
		lines(table$V1, table$V8)
		lines(table$V1, table$V9)
		lines(table$V1, table$V10)
		lines(table$V1, table$V11)
		lines(table$V1, table$V12)
		lines(table$V1, table$V13)
		lines(table$V1, table$V14)
		dev.off()

	This should have produced a plot named [`coal_rates_neooli.pdf`](res/coal_rates_neooli.pdf). This plot should look as shown below:<p align="center"><img src="img/coal_rates_neooli.png" alt="Relate" width="600"></p> As you can see, this plot shows far more variation, with some very young peaks and others that are even older than those seen with *Astatotilapia burtoni*. **Question 2:** What could have caused the strong peak in the coalescence rate around 6e+06, 6 million years ago? [(see answer)](#q1) Obviously, the format in which trees are saved in this file is not the Newick format. **Question 3:** Can you identify from the plot which peak shows the coalescence rates between *Neolamprologus olivaceous* and *Astatotilapia burtoni*?

In the next part of the tutorial, we will extract a set of the most reliable trees from the result files written by Relate. We will then generate a species trees from this set of "gene" trees, once again using ASTRAL.

* First, we will generate a subset of the results that contains those trees supported by a large number of mutations, for which we will set a minimum of 1,000 mutations. To generate this subset, use the following command:

		RelateExtract --mode RemoveTreesWithFewMutations --anc NC_031969.f5.sub1.phased.renamed.reestimated.anc --mut NC_031969.f5.sub1.phased.renamed.reestimated.mut --threshold 1000 -o NC_031969.f5.sub1.phased.renamed.reestimated.subset
	
* Have a look a the second line in the newly generated file `NC_031969.f5.sub1.phased.renamed.reestimated.subset.anc` to find out how many trees remain in the dataset after applying the above filter. To do so, you could use the following command:

		head -n 2 NC_031969.f5.sub1.phased.renamed.reestimated.subset.anc

* You should see that about 13,000 trees remain in the dataset. From these, we will randomly select 100 trees for conversion into Newick format. This could be done for a larger subset; however, the number of 100 trees will be sufficient for this tutorial. To convert 100 randomly selected trees into Newick format, use the following command:

		for i in `cat NC_031969.f5.sub1.phased.renamed.reestimated.subset.mut | tail -n +2 | cut -d ";" -f 2 | sort -R | head -n 100`
		do
			RelateExtract \
				--mode TreeAtSNPAsNewick \
				--anc NC_031969.f5.sub1.phased.renamed.reestimated.subset.anc \
				--mut NC_031969.f5.sub1.phased.renamed.reestimated.subset.mut \
				--bp_of_interest ${i} -o tree
		done

* Comine all trees in Newick format that were written in the previous step, using the following command:
	
		cat tree_at_*.newick > extracted.trees
	
* To clean up the directory, you could then remove the original tree files in Newick format, using this command:

		rm tree_at_*.newick

* Each of the trees in Newick has four tips per species, for the two haplotypes of each of the two individuals per species. Tips are labeled by numbers, which correspond to the order of species in the file [`NC_031969.f5.sub1.phased.renamed.sample`](res/NC_031969.f5.sub1.phased.renamed.sample) that was generated by RelateFileFormats earlier. To tell ASTRAL how haplotypes map to species, write the following text to a file named `samples.txt`:
	
		astbur: 0,1,2,3
		altfas: 4,5,6,7
		telvit: 8,9,10,11
		neobri: 12,13,14,15
		neocan: 16,17,18,19
		neochi: 20,21,22,23
		neocra: 24,25,26,27
		neogra: 28,29,30,31
		neohel: 32,33,34,35
		neomar: 36,37,38,39
		neooli: 40,41,42,43
		neopul: 44,45,46,47
		neosav: 48,49,50,51
		neowal: 52,53,54,55

* Finally, provide use the set of "gene" trees in file `extracted.trees` together with the list haplotypes and sample IDs in file `samples.txt` to ASTRAL to estimate the species tree. To do so, use the following command:
	
		java -jar astral.jar -i extracted.trees -a samples.txt -o species.tre

	This command should have resulted in the species tree shown below:<p align="center"><img src="img/figtree2.png" alt="FigTree" width="600"></p> This tree is largely congruent with the species trees determined with SVDQuartets in tutorial [Maximum-Likelihood Species Tree Inference](../ml_species_tree_inference/README.md) and SNAPP in tutorial [Divergence Time Estimation with SNP Data](../divergence_time_estimation_with_snp_data/README.md).

<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** Out dataset included 2 samples of each of 14 species, therefore 28 samples. For each sample, two phased sequences (= haplotypes) are included in the dataset, therefore there are 56 haplotypes in total.


<a name="q2"></a>

* **Question 2:** Given that also the first plot for *Astatotilapia burtoni* showed more variation in the rate estimates for older species that seemed to result only from stochasticity in the estimation, I suspect that the same could be responsible for the peak around 6 million years ago. Without further analyses, however, we can not exclude that it actually shows a true biological signal, such as possible introgression from a species that diverged as early as 6 million years ago from the other *Neolamprologus* species.


<a name="q3"></a>

* **Question 3:** The line for the coalescence rate between *Neolamprologus olivaceous* and *Astatotilapia burtoni* should be the only one that is present in both plots; therefore, it has to be the one that peaks just a bit later than 2e+06, 2 million years ago.