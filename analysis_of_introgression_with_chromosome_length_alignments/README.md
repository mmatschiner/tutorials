# Analysis of Introgression with Chromosome-Length Alignments

A tutorial on the analysis of hybridization and introgression with whole-chromosome alignments

## Summary

XXX

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Inferring recombination breakpoints with Saguaro](#saguaro)
* [Identifying alignment blocks for phylogenetic analysis](#alignments)
* [Automating BEAST2 analyses](#beast2)
* [Analyzing introgression with divergence times](#divtimes)
* [Simulating recombination with c-genie](#cgenie)


<a name="outline"></a>
## Outline

In this tutorial I am going XXX


<a name="dataset"></a>
## Dataset

XXX In contrast to the dataset used in tutorial [Analysis of Introgression with SNP Data](../analysis_of_introgression_with_snp_data/README.md), only a single individual per species will be used here to computational requirement of the phylogenetic analyses. XXX


<center>

| Sample ID | Species ID | Species name                  | Tribe         |
|-----------|------------|-------------------------------|---------------|
| IZC5      | astbur     | *Astatotilapia burtoni*       | Haplochromini |
| AUE7      | altfas     | *Altolamprologus fasciatus*   | Lamprologini  |
| JBD6      | telvit     | *Telmatochromis vittatus*     | Lamprologini  |
| JUH9      | neobri     | *Neolamprologus brichardi*    | Lamprologini  |
| LJC9      | neocan     | *Neolamprologus cancellatus*  | Lamprologini  |
| KHA7      | neochi     | *Neolamprologus chitamwebwai* | Lamprologini  |
| IVE8      | neocra     | *Neolamprologus crassus*      | Lamprologini  |
| JWH2      | neogra     | *Neolamprologus gracilis*     | Lamprologini  |
| JWG9      | neohel     | *Neolamprologus helianthus*   | Lamprologini  |
| JWH4      | neomar     | *Neolamprologus marunguensis* | Lamprologini  |
| JWH6      | neooli     | *Neolamprologus olivaceous*   | Lamprologini  |
| ISB3      | neopul     | *Neolamprologus pulcher*      | Lamprologini  |
| ISA8      | neosav     | *Neolamprologus savoryi*      | Lamprologini  |
| KFD2      | neowal     | *Neolamprologus walteri*      | Lamprologini  |

</center>


<a name="requirements"></a>
## Requirements

* **Saguaro:** [Saguaro](http://saguarogw.sourceforge.net) ([Zamani et al. 2013](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-347)) implements a hidden-Markov model for the detection of recombination breakpoints in chromosome-length alignments. Installation instructions can be found at [Saguaro's webpage](http://saguarogw.sourceforge.net), and the download is available from the associated [Sourceforge repository](https://sourceforge.net/p/saguarogw/code/HEAD/tree/) (click "Download Snapshot" to download). Note that Saguaro runs on Linux, and while in principle the installation should also be possible on a Mac, this does not seem to be easy and is not supported by the authors. The installation is not possible on Windows. If you have access to a Linux server with a Saguaro installation, but you would like to run the rest of the tutorial on your own machine, you can do so by transferring input and ouput files of Saguaro via scp between your machine and the Linux server. To check whether the Saguaro installation succeeded, just type `Saguaro` on the command line. If you should fail to install Saguaro, you can skip the tutorial part with the Saguaro analysis and continue with ready-made Saguaro result files afterwards.

* **RAxML:** If you followed tutorials [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md) and [Maximum-Likelihood Species-Tree Inference](../ml_species_tree_inference/README.md), you should have [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) ([Stamatakis 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)) installed on your machine already. If not, you will find source code for Mac OS X and Linux, as well as precompiled executables for Windows, on RAxML's github page [https://github.com/stamatak/standard-RAxML](https://github.com/stamatak/standard-RAxML). See tutorial [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md) for more details on the installation of RAxML.

* **BEAST2:** Like RAxML, you likely have the BEAST2 package installed already if you followed other tutorials in this collection. If not, you can downloaded it from the [BEAST2 website](https://www.beast2.org). As all these programs are written in Java, compilation is not required, and all programs should work on Mac OS X, Linux, and Windows. I recommend downloading the program versions that include the Java Runtime Environment, which may prevent conflicts with Java versions that may already be installed on your machine.<br>

* **Libraries for Python 3.x:** The Python libraries [dendropy](https://www.dendropy.org) ([Sukumaran and Holder 2010](https://academic.oup.com/bioinformatics/article/26/12/1569/287181)) and [msprime](https://msprime.readthedocs.io/en/stable/index.html) [(Kelleher et al. 2016)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004842) will be required for analyses of topology frequencies and for simulations of recombination, respectively. The two libraries can be installed with pip for Python 3, using the following commands:

		python3 -m pip install --user msprime
		python3 -m pip install --user dendropy
		
	The installations can be tested with these commands:
	
		python3 -c 'import msprime'
		python3 -c 'import dendropy'


* **Libraries for R:** Two R libraries will be required in this tutorial; these are [ape](https://cran.r-project.org/web/packages/ape/index.html) ([Paradis et al. 2004](https://academic.oup.com/bioinformatics/article/20/2/289/204981)) and [coda](https://cran.r-project.org/web/packages/coda/index.html) ([Plummer et al. 2006](http://oro.open.ac.uk/22547/)). If you're already familiar with R, just install these packages in the usual way. If not, the easiest way to do so might be via the command line. Type `R` to open the R environment interactively. Then, run the following commands:

		install.packages("ape", repos="http://ftp.gwdg.de/pub/misc/cran/", dependencies=T)
		install.packages("coda", repos="http://ftp.gwdg.de/pub/misc/cran/", dependencies=T)
		
	To ensure that both packages have been successfully installed, type these commands:
	
		library(ape)
		library(coda)

	If both commands result in no error messages, the two packages are ready to be used. Then, quit the R environment with `quit(save="no")`.



<a name="saguaro"></a>
## Inferring recombination breakpoints with Saguaro

In this part of the tutorial, the software [Saguaro](http://saguarogw.sourceforge.net) ([Zamani et al. 2013](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-347)) will be used to detect boundaries between chromosome regions that are characterized by different phylogenetic histories. However, for computational reasons, Saguaro does not infer these phylogenetic histories directly. Instead, the analysis performed by Saguaro is based on what the authors call "cacti", sets of distance matrices that describe how different each pair of genomes is relative to all others. These cacti can to some extent be considered as proxies for phylogenetic histories, as the genetic difference between pairs of chromosomes is obviously linked to their phylogenetic relatedness. Assuming (naively) that Saguaro accurately detects all recombination breakpoints, we applied this method in [Gante et al. (2016)](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13767) to select chromosomal blocks of 500 kb that were not broken up by recombination, and we used these blocks for subsequent phylogenetic analysis.

However, it has been shown based on simulations ([Martin and van Belleghem 2017](http://www.genetics.org/content/206/1/429.long)) that Saguaro may not always be able to discriminate chromosomal regions characterized by different phylogenetic histories and that subtle changes caused by incomplete lineage sorting may remain undetected. If this was the case for the chromosomal blocks used in [Gante et al. (2016)](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13767) for phylogenetic analysis, each of the inferred trees may have been affected by the issues that potentially result from concatenation, including exaggerated node support ([Kubatko and Degnan 2007](https://academic.oup.com/sysbio/article/56/1/17/1658327)) and overestimated divergence times ([Ogilvie et al. 2017](https://academic.oup.com/mbe/article/34/8/2101/3738283)).

To avoid recombination in the alignments used for phylogenetic analysis in this tutorial, we are therefore going to use Saguaro only as a first attempt of identifying recombination breakpoints. The breakpoints identied by Saguaro will be taken into account when generating alignment blocks for phylogenetic analysis; however, each of these alignment blocks will further be filtered according to a second test for recombination with the software [Phy Test](http://www.maths.otago.ac.nz/~dbryant/software.html) ([Bruen et al. 2006](http://www.genetics.org/content/172/4/2665)).

At the beginning of the analysis, Saguaro will calculate a single cactus for the entire alignment, and a score will be calcuated for each variable alignment position, describing the fit between this site and the first cactus. Based on these scores, genomic regions with a poor fit to the current cactus are identied with the hidden Markov model implemented in Saguaro, and a new cactus is defined for these. This process is repeated multiple times, thus further partitioning the alignment into segments, and at the same time assigning a cactus out of a growing set of cacti to each segment. Details of this procedure are described in ([Zamani et al. 2013](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-347)).

As in tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md) and other tutorials, we are here going to focus on chromosome 5 of tilapia only to reduce the computational demand of all analyses.

* To generate a chromosome-length sequence alignment for the 14 cichlid species included in the dataset, we'll need several different input files:

	* 	The VCF file `NC_031969.f5.sub1.phased.masked.mod.vcf.gz` with phased SNP variation generated in tutorial [Analysis of Introgression with SNP Data](../analysis_of_introgression_with_snp_data/README.md). Either copy this file from the analysis directory of that tutorial if you followed it already, or download file [`NC_031969.f5.sub1.phased.masked.mod.vcf.gz`](NC_031969.f5.sub1.phased.masked.mod.vcf.gz) by clicking on the link and copy it to the analysis directory for this tutorial.

	* As the VCF file `NC_031969.f5.sub1.phased.masked.mod.vcf.gz` does not contain information about invariant sites but these are required for phylogenetic inference with BEAST2, we assume that the sites between those included in the VCF file are in fact unchanged between the Lamprologini species and the tilapia reference sequence. Thus, we'll use the reference sequence for chromosome 5 of tilapia to fill the gaps between SNPs included in the VCF. To download the sequence for chromosome 5 of tilapia from GenBank, you can use the following command:

			wget 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&val=NC_031969&extrafeat=0&maxplex=1' -O NC_031969.fasta
		
		The downloaded sequence will in in Fasta format in a file named `NC_031969.fasta`.
		
	* The assumption that sites between SNPs are invariant is only justified for sites in which SNPs could have been called at all if there had been any. This is not the case for all sites between the SNPs, due to low sequencing coverage in some regions or low mapping quality in repetitive regions. It would also have been impossible to call SNPs very close to indel variation, because these were filtered when the original VCF file `NC_031969.f5.sub1.vcf.gz` (see tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md)) was prepared. Thus, we will take into account information about the chromosomal regions in which SNP variation could not have been called in the first place, and we will conservatively set these chromosomal regions to missing, coded by "N", in the generated alignment. The information about callable chromosomal regions is contained in a set of files in [BED format](http://genome.ucsc.edu/FAQ/FAQformat#format1); these files are within a compressed directory named [`masks.tgz`](data/masks.tgz). Download this directory by clicking on the link and copy it to your analysis directory.

* Unless the directory with the files in BED format has already been uncompressed automatically, do so with the following command:

		tar -xzf masks.tgz
		
* Have a look at the content of one of the files in BED format, for example using the following command:

		less masks/AUE7.NC_031969.f5.merged.bed
		
	You'll see that these files contain three columns, as shown in the text below:
	
		NC_031969       0       45144
		NC_031969       45169   88297
		NC_031969       88299   88363
		NC_031969       88365   88370
		NC_031969       88372   98367
		NC_031969       98368   98403
		NC_031969       98405   98409
		NC_031969       98411   98429
		NC_031969       98431   100023
		...
		
	The second and third of these columns indicate the beginning and the end of chromosomal regions that were masked because SNPs in these regions would not have been detected due to low coverage or other reasons. You'll see that these uncallable regions are quite long compared to the callable regions in between them. In the example shown above, only few callable sites were found within the first 100 kb of the chromosome. Overall, about two thirds of the chromosome are uncallable, which can be explained by a combination of several causes, including mapping to a distant reference, masking of repetitive regions, and rigorous filtering of sites close to insertions and deletions. However, even after conservatively setting all these regions to missing, the remaining sequence information will be more than sufficient for phylogenetic analysis.

* The three different types of information (the SNP data, the reference sequence, and the information on uncallable regions) will be combined to generate chromosome-length sequences in Fasta format by the Ruby script [`fill_seq.rb`](src/fill_seq.rb). To reduce the memory requirement of this script, we are going to run it separately for each of the samples included in the dataset, and as input for the script we are going to use individual VCF files that contain SNP data only for the respective sample. Prepare these individual VCF files from the full VCF file [`NC_031969.f5.sub1.phased.masked.mod.vcf.gz`](data/NC_031969.f5.sub1.phased.masked.mod.vcf.gz), using the following command (if you copied the VCF file from the analysis directory of tutorial [Analysis of Introgression with SNP Data](../analysis_of_introgression_with_snp_data/README.md), it might be uncompressed; if so, just remove the file extension `.gz` in this command):

		for i in IZC5 AUE7 JBD6 JUH9 LJC9 KHA7 IVE8 JWH2 JWG9 JWH4 JWH6 ISB3 ISA8 KFD2
		do
			bcftools view -s ${i} -o ${i}.NC_031969.f5.masked.mod.vcf NC_031969.f5.sub1.phased.masked.mod.vcf.gz
		done
		
	The above command should have written 14 separate uncompressed VCF files. To make sure that these have been generated, you could use `ls *.NC_031969.f5.masked.mod.vcf`.

* We can now run the script [`fill_seq.rb`](src/fill_seq.rb) to combine, for each sample, the SNP data with the reference sequence and the information on callable regions, using the following command (this will take about ten minutes to finish):

		for i in *.NC_031969.f5.masked.mod.vcf
		do
			sample_id=`basename ${i%.NC_031969.f5.masked.mod.vcf}`
			echo -n "Translating file ${i}..."
			ruby fill_seq.rb ${i} NC_031969.fasta masks/${sample_id}.NC_031969.f5.merged.bed ${sample_id}.NC_031969.f5.masked.mod.fasta
			echo " done."
		done

	As you can see from the command above, the script `fill_seq.rb` expects four command-line arguments. These are, in this order,
	
	* the VCF file with SNP information,
	* the reference sequence in Fasta format,
	* the file in BED format containing information about callable regions,
	* the name of an output file, which will be written in Fasta format.

	The script expects that all input files are for a single chromosome only, so we would have to run it multiple times if we would use information from several chromosomes.

* Have a look at one of the Fasta files generated by script `fill_seq.rb`, for example using the following command:

		less -S AUE7.NC_031969.f5.masked.mod.fasta
		
	You'll see that it contains two sequences, one for each haplotype of the phased VCF file.

* Next, combine all sequences into a single file in Fasta format named `NC_031969.f5.masked.mod.fasta`. Because all sequences are already aligned in the same way to the tilapia reference sequence, the combined file will already be perfectly aligned and no realignment with a tool like MAFFT (see tutorial [Multiple Sequence Alignment](../multiple_sequence_alignment/README.md)) will be necessary. The combined Fasta file `NC_031969.f5.masked.mod.fasta` will have a file size of about 1 Gb, thus make sure that you have sufficient disk space left before executing the following command:

		cat *.NC_031969.f5.masked.mod.fasta > NC_031969.f5.masked.mod.fasta
		
* You could now remove the per-sample Fasta files to save space, using the following command:

		rm *.NC_031969.f5.masked.mod.fasta

* As you've already seen, the Fasta file `NC_031969.f5.masked.mod.fasta` contains two sequences for each sample (and, given that we're here using a single sample per species, this means that it also contains two sequences per species). If we would use this dataset for analyses with Saguaro, these analyses might infer recombination breakpoints at which only the relationships between the two sequences per species change but not the relationships among species. While it would be a more conservative approach to also account for these recombination breakpoints, one could argue that their influence on the phylogenetic inference of the species tree is likely minimal. Therefore, and in order to shorten the run time of the Saguaro analysis, we are going to reduce the dataset for the Saguaro analysis to a single sequence per species, meaning that only recombination breakpoints affecting the relationships among species will be detected. To write the first sequence of each sample to a new file named `NC_031969.f5.masked.mod.A.fasta`, use the following command:

		cat NC_031969.f5.masked.mod.fasta | grep -A 1 "_A" | grep -v "^--$" > NC_031969.f5.masked.mod.A.fasta

* The next command can only be executed if Saguaro has been installed. Thus, if you could not install Saguaro on your own computer but on a remote Linux server, you should now copy the file `NC_031969.f5.masked.mod.A.fasta` to the Linux server, with the `scp` or `rsync` command-line tools.

* To prepare the sequence alignment for analysis with Saguaro, it needs to be converted from Fasta format into the "Feature" format of Saguaro. This can be done with the tool Fasta2HMMFeature that is part of the Saguaro installation. To see the available options of Fasta2HMMFeature, just type the program name on the command line:

		Fasta2HMMFeature
		
	This will show the following text:
	
		Fasta2HMMFeature: Converts multi-fasta data into a Saguaro-digestable file of features

		Available arguments:
		
		-i<string> : input fasta file (multiple alignment)
		-o<string> : binary output file
		-nosame<bool> : skip positions in which all calls are the same (def=0)
		-m<int> : minimum coverage (def=2)
		-d<int> : minimum disagreeing (def=2)
		-n<string> : chromosome name (def=mult)

	With the options `-m` and `-d` we can specify that Saguaro should only use positions with a minimum coverage (here, this means non-missing sequence information) or to use only those where at least two sequences are different from the others. The default values for both is 2, which means that only parsimony-informative sites are considered.
	
* Use the Fasta file `NC_031969.f5.masked.mod.A.fasta` as the input for the Saguaro analysis with option `-i` and name the output `NC_031969.f5.masked.mod.A.feature` with option `-o`. Also make sure to set the chromosome name to `NC_031969.f5.masked.mod` (without the `.A` at the end) with option `-n`; this will be important when the Saguaro results are later used to generate alignment blocks. Thus, execute Fasta2HMMFeature with the following command:

		Fasta2HMMFeature -i NC_031969.f5.masked.mod.A.fasta -o NC_031969.f5.masked.mod.A.feature -n NC_031969.f5.masked.mod

* Next, have a look at the available options for Saguaro, by typing the program name on the command line:

		Saguaro
		
	This should output the following text:
	
		Saguaro: Smoothly, automatically and generically uncover the ancestry of related organisms.


		Available arguments:

		-f<string> : Feature vector (def=)
		-l<string> : Feature vector list file (def=)
		-o<string> : output directory
		-cycle<int> : iterations per cycle (def=2)
		-iter<int> : iterations with split (def=40)
		-t<double> : transition penalty (def=150)
		-resume<int> : resume w/ iteration # (def=0)
		-neurons<int> : number of neurons in the SOM (def=800)

		
	You'll see that input can either be provided with option `-f` or option `-l`. Of these, `-f` followed by the name of a Feature file should be used when you have just a single such file (as in our case), and a list of file names can be provided with with `-l` if you have multiple Feature files for different chromosomes. The options `-cycle`, `-iter`, and `-neurons` determine settings for the hidden Markov model of Saguaro. The default values for these settings are relatively thorough and usually appropriate. To reduce the run time of Saguaro, we will here, however, use less cycles and iterations as well as a smaller number of neurons. Another parameter that might be worth modifying is the transition penalty, specified with option `-t`. By reducing the transition penalty, the sensitivity of Saguaro to detect recombination breakpoints can be increased. This means that Saguaro would then be less likely to miss true recombination breakpoints but it might also lead to the detection of false positives. To run Saguaro with `NC_031969.f5.masked.mod.A.feature` as the input file, `saguaro_results` as the name of the output directory, 10 iterations with 1 cycle per iterations, and 100 neurons, use the following command:
	
		Saguaro -f NC_031969.f5.masked.mod.A.feature -o saguaro_results -iter 10 -cycle 1 -neurons 100

* Once the Saguaro analysis has finished, results should have been written to directory `saguaro_results`. Have a look at the content of this directory. You'll see the following files:

		HMMTrain.out.0
		HMMTrain.out.1
		HMMTrain.out.2
		...
		LocalTrees.out
		saguaro.cactus
		saguaro.config
		saguaro.garbage
		saguaro.garbage.vec
		
	Of these, only `LocalTrees.out` and `saguaro.cactus` are relevant for us. If you ran the Saguaro analysis on a Linux server because you could not install Saguaro on your own computer, you could now copy these files back to your own computer with `scp` or `rsync`. If you were unable to run the Saguaro analysis at all, you can download the ready-made result files [`LocalTrees.out`](res/LocalTrees.out) and [`saguaro.cactus`](res/saguaro.cactus) with the two links.

* Open file `saguaro.cactus` in a text editor, or from the command line with the following command:

		less saguaro.cactus

	You should see more or less the following content:
	
		cactus0
		AUE7_A  ISA8_A  ISB3_A  IVE8_A  IZC5_A  JBD6_A  JUH9_A  JWG9_A  JWH2_A  JWH4_A  JWH6_A  KFD2_A  KHA7_A  LJC9_A
		AUE7_A  -0.000000       0.793221        0.872338        0.839333        0.524633        0.687894        0.843126        0.866756        0.819825        0.839036        0.877009        0.834673        0.834047        0.069426
		ISA8_A  0.793221        -0.000000       0.304833        0.368500        0.444325        0.384680        0.335022        0.323023        0.322399        0.358284        0.331484        0.486985        0.495132        0.906134
		ISB3_A  0.872338        0.304833        -0.000000       0.356743        0.493709        0.430859        0.232075        0.288845        0.309689        0.357865        0.222982        0.526403        0.524468        0.987994
		IVE8_A  0.839333        0.368500        0.356743        -0.000000       0.478407        0.410530        0.337957        0.340104        0.314575        0.189210        0.340049        0.490898        0.504986        0.935056
		IZC5_A  0.524633        0.444325        0.493709        0.478407        -0.000000       0.390289        0.503307        0.485551        0.459490        0.475014        0.496944        0.486223        0.500198        0.650715
		JBD6_A  0.687894        0.384680        0.430859        0.410530        0.390289        -0.000000       0.427936        0.410087        0.411593        0.412935        0.439920        0.422403        0.432646        0.782539
		JUH9_A  0.843126        0.335022        0.232075        0.337957        0.503307        0.427936        -0.000000       0.261218        0.301458        0.370035        0.287994        0.537055        0.533601        0.951846
		JWG9_A  0.866756        0.323023        0.288845        0.340104        0.485551        0.410087        0.261218        -0.000000       0.299211        0.368540        0.292421        0.512795        0.524726        0.958697
		JWH2_A  0.819825        0.322399        0.309689        0.314575        0.459490        0.411593        0.301458        0.299211        -0.000000       0.313155        0.343014        0.504644        0.509835        0.921273
		JWH4_A  0.839036        0.358284        0.357865        0.189210        0.475014        0.412935        0.370035        0.368540        0.313155        -0.000000       0.360356        0.475263        0.494985        0.951263
		JWH6_A  0.877009        0.331484        0.222982        0.340049        0.496944        0.439920        0.287994        0.292421        0.343014        0.360356        -0.000000       0.525858        0.537397        0.963240
		KFD2_A  0.834673        0.486985        0.526403        0.490898        0.486223        0.422403        0.537055        0.512795        0.504644        0.475263        0.525858        -0.000000       0.079350        0.938385
		KHA7_A  0.834047        0.495132        0.524468        0.504986        0.500198        0.432646        0.533601        0.524726        0.509835        0.494985        0.537397        0.079350        -0.000000       0.929866
		LJC9_A  0.069426        0.906134        0.987994        0.935056        0.650715        0.782539        0.951846        0.958697        0.921273        0.951263        0.963240        0.938385        0.929866        -0.000000
		cactus1
		...

	This is the distance matrix for the first cactus (cactus0), followed by distance matrices for all other inferred cacti.

* Next, open file [`LocalTrees.out`](res/LocalTrees.out) in a text editor or on the command line with the `less` command. You'll see something like this:

		Reading features...
		done!
		Reading models.
		Dynprog'ing...
		Setting up HMM...
		Adding word cactus0 as # 0
		...

	Scroll down a bit to this part:

		...
		Processed: 80000 (95.6034 %)
		REPORTING Traceback and Update
		cactus7 NC_031969.f5.masked.mod: 204750 - 3479153       length: 3274403 (frames 1-2704 l=2703)  score=62.9241
		AUE7_A  ISA8_A  ISB3_A  IVE8_A  IZC5_A  JBD6_A  JUH9_A  JWG9_A  JWH2_A  JWH4_A  JWH6_A  KFD2_A  KHA7_A  LJC9_A
		AUE7_A  -0.00   0.77    0.92    0.88    0.65    0.69    0.89    0.87    0.87    0.84    0.88    0.80    0.76    0.09
		ISA8_A  0.77    -0.00   0.46    0.47    0.59    0.49    0.42    0.42    0.43    0.45    0.42    0.53    0.55    0.87
		ISB3_A  0.92    0.46    -0.00   0.34    0.55    0.52    0.30    0.40    0.34    0.40    0.24    0.52    0.55    0.97
		IVE8_A  0.88    0.47    0.34    -0.00   0.57    0.49    0.34    0.39    0.36    0.25    0.38    0.56    0.57    0.92
		IZC5_A  0.65    0.59    0.55    0.57    -0.00   0.53    0.61    0.63    0.61    0.61    0.62    0.54    0.54    0.74
		JBD6_A  0.69    0.49    0.52    0.49    0.53    -0.00   0.52    0.51    0.53    0.52    0.51    0.46    0.48    0.73
		JUH9_A  0.89    0.42    0.30    0.34    0.61    0.52    -0.00   0.38    0.33    0.42    0.33    0.53    0.58    0.95
		JWG9_A  0.87    0.42    0.40    0.39    0.63    0.51    0.38    -0.00   0.37    0.39    0.39    0.54    0.57    0.94
		JWH2_A  0.87    0.43    0.34    0.36    0.61    0.53    0.33    0.37    -0.00   0.42    0.37    0.57    0.56    0.95
		JWH4_A  0.84    0.45    0.40    0.25    0.61    0.52    0.42    0.39    0.42    -0.00   0.43    0.53    0.55    0.91
		JWH6_A  0.88    0.42    0.24    0.38    0.62    0.51    0.33    0.39    0.37    0.43    -0.00   0.55    0.56    0.94
		KFD2_A  0.80    0.53    0.52    0.56    0.54    0.46    0.53    0.54    0.57    0.53    0.55    -0.00   0.16    0.85
		KHA7_A  0.76    0.55    0.55    0.57    0.54    0.48    0.58    0.57    0.56    0.55    0.56    0.16    -0.00   0.84
		LJC9_A  0.09    0.87    0.97    0.92    0.74    0.73    0.95    0.94    0.95    0.91    0.94    0.85    0.84    -0.00
		cactus4 NC_031969.f5.masked.mod: 3481523 - 3581323      length: 99800   (frames 2705-2808 l=103)        score=1692.5
		...
		
	These lines contain information on the first segment identified by Saguaro. It is assigned to cactus7, and is located between position 204750 and position 3479153 on the chromosome named "NC_031969.f5.masked.mod". In contrast to distance matrices given in file `saguaro.cactus`, the distance matrices given in this in `LocalTrees.out` are calculated per segment, not per cactus. Nevertheless, the distance matrix of the first segment between position 204750 and position 3479153 is likely very similar to that of cactus7, otherwise, this segment would not have been assigned to this cactus.

	As we're not particularly interested in the distance matrices of segments, but more in the placement of segment boundaries so that we can select alignment regions for phylogenetic analyses that are not broken up by boundaries, the most imporant information for us is in the header lines for each segment.

* To see only header information for each cactus, type this command:

		cat LocalTrees.out | grep length
		
	You should see output like this:

		cactus7	NC_031969.f5.masked.mod: 204750 - 3479153	length: 3274403	(frames 1-2704 l=2703) 	score=62.9241
		cactus4	NC_031969.f5.masked.mod: 3481523 - 3581323	length: 99800	(frames 2705-2808 l=103) 	score=1692.5
		cactus7	NC_031969.f5.masked.mod: 3590084 - 6158718	length: 2568634	(frames 2809-8397 l=5588) 	score=104.99
		cactus4	NC_031969.f5.masked.mod: 6158744 - 6166115	length: 7371	(frames 8398-8429 l=31) 	score=18425
		cactus5	NC_031969.f5.masked.mod: 6166197 - 6167838	length: 1641	(frames 8430-8438 l=8) 	score=65580.2
		cactus6	NC_031969.f5.masked.mod: 6167956 - 6182408	length: 14452	(frames 8439-8507 l=68) 	score=8637.73
		cactus5	NC_031969.f5.masked.mod: 6182675 - 6188272	length: 5597	(frames 8508-8523 l=15) 	score=37313
		...

* To find out how many segments Saguaro has identified, you could type the following command:

		cat LocalTrees.out | grep length | wc -l

* In order to visualize recombination breakpoints and the cacti assigned to the different chromosomal regions, you can use the Ruby script [`paint_chromosomes.rb`](src/paint_chromosomes.rb). To find out how to run this script, you could type this command:

		ruby paint_chromosomes.rb

	This should show the following help text:

		paint_chromosomes.rb

		This script uses output from the software Saguaro to paint
		chromosomes according to Saguaro cacti. The output will be
		written in svg format.

		This script should be run e.g. with
		ruby paint_chromosomes.rb LocalTrees.out LocalTrees.svg
		where 'LocalTrees.out' should be replaced with the actual path
		to the Saguaro output file.

* Now, run the script `paint_chromosomes.rb` with the following command:

		ruby paint_chromosomes.rb LocalTrees.out
		
	This will generate a vector graphic file in SVG format that will be written to the same directory in which `LocalTrees.out` is placed, and it will be named `LocalTrees.svg`.
	
* Open the vector graphic file [`LocalTrees.svg`](res/LocalTrees.svg) with a program capable of reading files in SVG format, for example with a browser such as Firefox or with Adobe Illustrator. This plot should look as shown below:<p align="center"><img src="img/LocalTrees.png" alt="Saguaro" width="800"></p>

	In this visualization, chromosomal regions assigned to the most common cactus are drawn in dark gray, and segments assigned to other cacti are shown in red, orange, cyan, and light green, purple (in decreasing frequency). With more than six different cacti, all remaining cacti are shown in light gray. You'll notice that the most frequent cactus (in dark gray) is found in the first part of the chromosome, that the last part of the chromosome is dominated by the third-most frequent cactus (in orange), and that the regions in the center of the chromosome are generally shorter than those at the chromosome ends. This observation disagrees with the pattern expected if recombination was more frequent in in the peripheries of the chromosome, which is commonly observed in animals and plants ([Haenel et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14699)). But as we will see below, the alignment contains a higher amount of missing data towards the chromosome ends, which is likely responsible for the lower number of detected recombination breakpoints those parts of the chromosome. However, keep in mind that these results may not be particularly accurate, as we've only performed a short Saguaro analysis.

<a name="alignments"></a>
## Identifying alignment blocks for phylogenetic analysis

After having inferred recombination breakpoints with Saguaro, we can now cut the breakpoint-free parts of the chromosome-length alignment into alignment blocks of a fixed length, and we can try to identify the most suitable of these blocks for subsequent phylogenetic analysis. Ideally, the blocks used for phylogenetic analysis should have as little missing data as possible, be as informative as possible, and show no signs of recombination that may have remained undetected in the Saguaro analysis. In this part of the tutorial, we are going to quantify these characteristics for all alignment blocks to select a set of blocks for phylogenetic analysis with BEAST2 in the next part of the tutorial.

As the length of each block, we here use 50 kb, assuming that this length is a good compromise between increasing probability of undetected recombination with longer blocks and decreasing phylogenetic signal with shorter blocks. For a more thorough analysis, however, it might be worth testing this assumption with blocks of different sizes.

* To cut the breakpoint-free parts of the alignment into individual alignment blocks, we can use the Ruby script `generate_alignments.rb`. This script was specifically written to read output by Saguaro as well as the input file(s) used in the Saguaro analysis. It expects four command-line arguments; these are
	* the name of the `LocalTrees.out` output file of Saguaro,
	* the name of a directory in which the alignment files used for the Saguaro analysis are located,
	* the name of a new directory to which all output files will be written,
	* the length of the alignment blocks.

	Run the script with the following command:
	
		ruby generate_alignments.rb LocalTrees.out . alignment_blocks 50000
		
	This will generate as many non-overlapping alignment blocks of 50 kb as possible without including the recombination breakpoints identified by Saguaro. With the settings used in the command above, the script will read file `LocalTrees.out` and will identify chromosome names based on the information in this file. In our case, it will identify "NC_031969.f5.masked.mod" as a chromosome name, and it will search for file `NC_031969.f5.masked.mod.fasta` in the directory that's specified as the second command line argument, which here simply is `.`, the shortcut for the current directory. Output files in Nexus format will be written to a directory named `alignment_blocks`, which will be created inside the current directory.

* Have a look at the new directory named `alignment_blocks`. Note that the files in this directory are named according to the name of the linkage group and the first and the last position of the alignment block:

		NC_031969.f5.masked.mod_00666951_00716950.nex
		NC_031969.f5.masked.mod_01316951_01366950.nex
		NC_031969.f5.masked.mod_01366951_01416950.nex
		NC_031969.f5.masked.mod_01566951_01616950.nex
		NC_031969.f5.masked.mod_01716951_01766950.nex
		NC_031969.f5.masked.mod_01816951_01866950.nex
		NC_031969.f5.masked.mod_02066951_02116950.nex
		...

* Find out how many alignment blocks were written, using the following command:

		ls alignment_blocks/*.nex | wc -l
		
* To facilitate the phylogenetic analyses, it will help to simplify the alignment names with the following set of commands.

		for i in alignment_blocks/*.nex
		do
			new_name=`echo ${i} | sed 's/\.f5\.masked\.mod//g'`
			mv ${i} ${new_name}
		done

* The files in directory `alignment_blocks` should now be named like this:

		NC_031969_00666951_00716950.nex
		NC_031969_01316951_01366950.nex
		NC_031969_01366951_01416950.nex
		NC_031969_01566951_01616950.nex
		NC_031969_01716951_01766950.nex
		NC_031969_01816951_01866950.nex
		NC_031969_02066951_02116950.nex
		...

* Pick one of the files in this directory at random, and open it either in a text editor or in AliView just to get a feeling for the size of the alignment block, as well as for its sequence variation and the amount of missing data. You should see something like this:<p align="center"><img src="img/aliview1.png" alt="AliView" width="600"></p>

	Like the alignment shown above, many alignments contain a very large proportion of missing data due to the strict filtering based on coverage, mapping quality, and indel proximity. Some alignment may even consist exclusively of missing data.

* We are going to use mean overall node support resulting from maximum-likelihood tree inference as one measure of phylogenetic informativeness. This inference will be done with the software RAxML, which you will be familiar with if you followed the tutorials [Maximum-Likelihood Phylogenetic Inference](../ml_phylogeny_inference/README.md) or [Maximum-Likelihood Species-Tree Inference](../ml_species_tree_inference/README.md). If you didn't, you may want to look up the instructions for these tutorials to learn more about RAxML analyses. However, to avoid RAxML analyses of alignment blocks without sequences data or with very little sequence data, we are first going to remove those alignment blocks with a proportion of missing data above 90%. We can do so with the following set of commands that uses the Ruby script [`get_proportion_of_missing_data.rb`](src/get_proportion_of_missing_data.rb). The input and output for this script is quite simple: It the name of a file in Nexus format as input and reports nothing but the proportion of missing data on the screen. With the set of commands below, this proportion is assigned to the variable "p\_missing" and if this variable is above 0.9, the Nexus-format alignment file is removed. Execute this set of commands:

		for i in alignment_blocks/*.nex
		do
			p_missing=`ruby get_proportion_of_missing_data.rb ${i}`
			if [[ ${p_missing} > 0.9 ]]
			then
				rm ${i}
			fi
		done

* Now, count the number of alignment files in directory `alignment_blocks` once again with the following command:

		ls alignment_blocks/*.nex | wc -l

	**Question 1:** How many alignments were removed due to a proportion of missing data greater than 90%? [(see answer)](#q1)

* Next, we are going to run maximum-likelihood tree inference with the software RAxML for each alignment file. We are going to do so with the following set of commands that first of all uses the Python script `convert.py` to convert each alignment file from Nexus to Phylip format, as RAxML does not accept Nexus format. The RAxML analysis will be performed with the "GTRCAT" model of [Stamatakis 2006](https://ieeexplore.ieee.org/abstract/document/1639535/); however, as were are now only interested in the mean node support and not in the inferred tree itself, the choice of substitution model is not going to have a strong influence on our results. By using the options `-f a -x ${RANDOM} -N 100`, we specify that the maximum-likelihood analysis should be followed by bootstrapping to assess node support with 100 replicates. With the commands below, RAxML is set to use four CPUs with the option `-T 4` to speed up the analysis, but if you have less CPUs available on your machine, you should change this number to e.g. `-T 2`. Then, run this set of commands:

		for i in alignment_blocks/*.nex
		do
			# Get the alignment id.
			id=`basename ${i%.nex}`

			# Convert the alignment to Phylip format.
			python3 convert.py -f phylip ${i} tmp.phy
			
			# Run maximum-likelihood tree inference and assess node support by bootstrapping.
			raxml -T 4 -m GTRCAT -n ${id} -s tmp.phy --silent -p ${RANDOM} -f a -x ${RANDOM} -N 100
			
			# Clean up RAxML output files.
			rm tmp.phy*
			rm RAxML_bipartitionsBranchLabels.*
			rm RAxML_bestTree.*
			rm RAxML_bootstrap.*
			mv RAxML_bipartitions.* alignment_blocks/${id}.tre
			mv RAxML_info.* alignment_blocks/${id}.info
		done
		
	With the about 470 alignment blocks and bootstrapping to be performed for each of them, these RAxML analyses will take a rather long time, around 7 hours. Thus, unless you have the chance to run these analyses overnight, you'll probably rather want to continue the rest of this tutorial with the ready-made RAxML results from my analysis. 

* Get summary statistics for each alignment.

		echo -e "block_id\tp_missing\tn_pi_sites\tmean_bootstrap_support\tphi_p" > block_stats.txt
		for i in alignment_blocks/*.nex
		do
			# Get the block ID.
			block_id=`basename ${i%.nex}`
			
			# Provide feedback on screen.
			echo -n "Analyzing file ${i}..."
			
			# Get the proportion of missing data.
			proportion_of_missing_data=`ruby get_proportion_of_missing_data.rb $i`

			# Get the number of parsimony-informative sites.
			number_of_pi_sites=`ruby get_number_of_pi_sites.rb $i`
						
			# Get the mean bootstrap support.
			bootstrap_value=`python3 get_mean_node_support.py ${i%.nex}.tre`
			bootstrap_value_strip=`echo -n $bootstrap_value`
			
			# Run PhiPack to test for recombination and get the resulting p-value.
			mkdir ${block_id}
			python3 convert.py ${i} -f fasta > ${block_id}/tmp.fasta
			cd ${block_id}
			p_value=`Phi -f tmp.fasta | tail -n 2 | head -n 1 | cut -d ":" -f 2 | tr -d '[[:space:]]'`
			cd ../
			rm -r ${block_id}
			
			# Print all output.
			echo -e -n ${block_id} >> block_stats.txt
			echo -e -n "\t" >> block_stats.txt
			printf "%.3f" ${proportion_of_missing_data} >> block_stats.txt
			echo -e -n "\t" >> block_stats.txt
			echo -n ${number_of_pi_sites} >> block_stats.txt
			echo -e -n "\t" >> block_stats.txt
			printf "%.1f" ${bootstrap_value_strip} >> block_stats.txt
			echo -e -n "\t" >> block_stats.txt
			echo -e -n ${p_value} >> block_stats.txt
			echo >> block_stats.txt
			
			# Provide feedback on screen.
			echo " done."

		done

* Make some plots in R to see correlations between the different parameters.

		R
	
	Then, ...
		
		t <- read.table("block_stats.txt", header=T)
		get_third_as_num <- function(x, split="_"){ return(as.numeric(strsplit(x, split=split)[[1]][3])) }
		get_fourth_as_num <- function(x, split="_"){ return(as.numeric(strsplit(x, split=split)[[1]][4])) }
		froms <- unlist(lapply(as.character(t$block_id), get_third_as_num))
		tos <- unlist(lapply(as.character(t$block_id), get_fourth_as_num))
		block_centers <- (froms + tos)/2
		
		pdf("alignment_blocks_missing.pdf", height=7, width=7)
		plot(block_centers, t$p_missing, xlab="Position", ylab="Proportion missing", main="NC_031969")
		abline(h=0.8)
		rect(-5000000, 0, 40000000, 0.8, col=rgb(0.164,0.629,0.594,alpha=0.25))
		dev.off()

	<p align="center"><img src="img/alignment_blocks_missing.png" alt="Number of PI sites" width="600"></p>

		pdf("alignment_blocks_pi_sites.pdf", height=7, width=7)
		plot(t$p_missing, t$n_pi_sites, xlab="Proportion missing", ylab="Number PI sites", main="NC_031969")
		abline(h=250)
		abline(v=0.8)
		rect(0, 250, 0.8, 1000, col=rgb(0.164,0.629,0.594,alpha=0.25))
		dev.off()
		
	<p align="center"><img src="img/alignment_blocks_pi_sites.png" alt="Number of PI sites" width="600"></p>
		
		pdf("alignment_blocks_bootstrap.pdf", height=7, width=7)
		plot(t$p_missing, t$mean_bootstrap_support, xlab="Proportion missing", ylab="Mean Bootstrap support", main="NC_031969")
		abline(h=70)
		abline(v=0.8)
		rect(0, 70, 0.8, 1000, col=rgb(0.164,0.629,0.594,alpha=0.25))
		dev.off()
		quit(save="no")

	<p align="center"><img src="img/alignment_blocks_bootstrap.png" alt="Mean bootstrap support" width="600"></p>


* Test some parameter combinations for filtering, how many blocks are left then? For example:

		ruby filter_blocks.rb block_stats.txt 0.8 250 70 0.05 | wc -l
		
* Make a list of filtered blocks:

		ruby filter_blocks.rb block_stats.txt 0.8 250 70 0.05 > block_stats_filtered.txt
			
* Copy the filtered blocks to a new directory.

		for i in `cat block_stats_filtered.txt | tail -n +2 | cut -f 1`
		do
			mkdir -p alignment_blocks_filtered/${i}
			cp alignment_blocks/${i}.nex alignment_blocks_filtered/${i}
		done


<a name="beast2"></a>
## Automating BEAST2 analyses

XXX

* Make a constraints file. Save the following as a new file named `constraints.xml`:

						<distribution id="All.prior" monophyletic="true" spec="beast.math.distributions.MRCAPrior" tree="@tree.t:Species">
							<taxonset id="All" spec="TaxonSet">
								<taxon idref="IZC5_A"/>
								<taxon idref="IZC5_B"/>
								<taxon idref="AUE7_A"/>
								<taxon idref="AUE7_B"/>
								<taxon idref="JBD6_A"/>
								<taxon idref="JBD6_B"/>
								<taxon idref="JUH9_A"/>
								<taxon idref="JUH9_B"/>
								<taxon idref="LJC9_A"/>
								<taxon idref="LJC9_B"/>
								<taxon idref="KHA7_A"/>
								<taxon idref="KHA7_B"/>
								<taxon idref="IVE8_A"/>
								<taxon idref="IVE8_B"/>
								<taxon idref="JWH2_A"/>
								<taxon idref="JWH2_B"/>
								<taxon idref="JWG9_A"/>
								<taxon idref="JWG9_B"/>
								<taxon idref="JWH4_A"/>
								<taxon idref="JWH4_B"/>
								<taxon idref="JWH6_A"/>
								<taxon idref="JWH6_B"/>
								<taxon idref="ISB3_A"/>
								<taxon idref="ISB3_B"/>
								<taxon idref="ISA8_A"/>
								<taxon idref="ISA8_B"/>
								<taxon idref="KFD2_A"/>
								<taxon idref="KFD2_B"/>
							</taxonset>
							<LogNormal meanInRealSpace="true" name="distr" offset="0.0">
								<parameter estimate="false" name="M">6.2666</parameter>
								<parameter estimate="false" name="S">0.13</parameter>
							</LogNormal>
						</distribution>
						
						<distribution id="Outgroup.prior" monophyletic="true" spec="beast.math.distributions.MRCAPrior" tree="@tree.t:Species">
							<taxonset id="Outgroup" spec="TaxonSet">
								<taxon idref="IZC5_A"/>
								<taxon idref="IZC5_B"/>
							</taxonset>
						</distribution>

						<distribution id="Ingroup.prior" monophyletic="true" spec="beast.math.distributions.MRCAPrior" tree="@tree.t:Species">
							<taxonset id="Ingroup" spec="TaxonSet">
								<taxon idref="AUE7_A"/>
								<taxon idref="AUE7_B"/>
								<taxon idref="JBD6_A"/>
								<taxon idref="JBD6_B"/>
								<taxon idref="JUH9_A"/>
								<taxon idref="JUH9_B"/>
								<taxon idref="LJC9_A"/>
								<taxon idref="LJC9_B"/>
								<taxon idref="KHA7_A"/>
								<taxon idref="KHA7_B"/>
								<taxon idref="IVE8_A"/>
								<taxon idref="IVE8_B"/>
								<taxon idref="JWH2_A"/>
								<taxon idref="JWH2_B"/>
								<taxon idref="JWG9_A"/>
								<taxon idref="JWG9_B"/>
								<taxon idref="JWH4_A"/>
								<taxon idref="JWH4_B"/>
								<taxon idref="JWH6_A"/>
								<taxon idref="JWH6_B"/>
								<taxon idref="ISB3_A"/>
								<taxon idref="ISB3_B"/>
								<taxon idref="ISA8_A"/>
								<taxon idref="ISA8_B"/>
								<taxon idref="KFD2_A"/>
								<taxon idref="KFD2_B"/>
							</taxonset>
						</distribution>


* Generate XML input files for BEAST2.

		for i in alignment_blocks_filtered/*
		do
			run_id=`basename ${i}`
			ruby beauti.rb -id ${run_id} -n ${i} -o ${i} -c constraints.xml -l 500000 -m HKY
		done
		
* Run all BEAST2 analyses.

		home=`pwd`
		for i in alignment_blocks_filtered/*
		do
			block_id=`basename ${i}`
			cd ${i}
			/Applications/Beast/2.5.0/bin/beast ${block_id}.xml
			cd ${home}
		done

* Open R:

		R
		
* Then:

		# Load the coda library.
		library(coda)
		
		# Read the log file as a table.
		t <- read.table("alignment_blocks_filtered/NC_031969_02066951_02116950/NC_031969_02066951_02116950.log", header=T)
		
		# Get the chain length as the number of rows in the table.
		chain_length = dim(t)[1]
		
		# Convert the table into a mcmc object.
		MCMCdata = mcmc(data=t, start=1, end=chain_length, thin=1)
		
		# Make a subset of the mcmc object, removing the burning and some parameter traces.
		MCMCsub <- as.mcmc(MCMCdata[(chain_length/10):chain_length, c(2, 3, 4, 8, 9, 10)])
		
		# Show a summary of the parameter estimates.
		summary(MCMCsub)
		
		# Calculate effective sample sizes for all parameters.
		effectiveSize(MCMCsub)
		
		# Find the lowest effective sample size.
		min(effectiveSize(MCMCsub))
		
* Then:

		# Make a pdf plot of the selected parameter traces.
		pdf("mcmc_traces.pdf", height=7, width=7)
		plot(MCMCsub, trace=TRUE, density=FALSE, smooth=TRUE, auto.layout=TRUE)
		dev.off()
		quit(save="no")

	<p align="center"><img src="img/mcmc_traces.png" alt="Coda" width="600"></p>

* Write the following code to a new file named `get_min_ess.r`:

		# Load the coda library.
		library(coda)
		
		# Get the command-line arguments.
		args <- commandArgs(trailingOnly = TRUE)
		log_file_name <- args[1]

		# Read the log file as a table.
		t <- read.table(log_file_name, header=T)
		
		# Get the chain length as the number of rows in the table.
		chain_length = dim(t)[1]
		
		# Convert the table into a mcmc object.
		MCMCdata = mcmc(data=t, start=1, end=chain_length, thin=1)

		# Make a subset of the mcmc object, removing the burning and some parameter traces.
		MCMCsub <- as.mcmc(MCMCdata[(chain_length/10):chain_length, c(2, 3, 4, 8, 9, 10)])
		
		# Find the lowest effective sample size.
		cat(min(effectiveSize(MCMCsub)), "\n")
		
* Then, run the R script first with a single file, the log file from the first alignment named [`NC_031969_02066951_02116950.log`](res/NC_031969_02066951_02116950.log):

		Rscript get_min_ess.r alignment_blocks_filtered/NC_031969_02066951_02116950/NC_031969_02066951_02116950.log

* Then, use the script to find the lowest ESS value of each log file:

		for i in alignment_blocks_filtered/*/*.log
		do
			echo -n -e "${i}\t"
			Rscript get_min_ess.r ${i} | tail -n 1
		done

* Remove the directories of all blocks with ESS values below 100:

		for i in alignment_blocks_filtered/*/*.log
		do
			block_id=`basename ${i%.log}`
			min_ess=`Rscript get_min_ess.r ${i} | tail -n 1 | cut -d "." -f 1`
			if [[ ${min_ess} -lt 100 ]]
			then
				rm -r alignment_blocks_filtered/${block_id}
				echo "Deleted directory alignment_blocks_filtered/${block_id}."
			fi
		done
		
* Make MCC summary trees for all genes:

		for i in alignment_blocks_filtered/*/*.trees
		do
			/Applications/Beast/2.5.0/bin/treeannotator -burnin 10 -heights mean ${i} ${i%.trees}.tre
		done

* Check the mean node support of all trees, using Python script [`get_mean_node_support.py`](src/get_mean_node_support.py):

		for i in alignment_blocks_filtered/*/*.tre
		do
			echo -ne "${i}\t"
			python3 get_mean_node_support.py ${i}
		done

* Combine all MCC trees in a single file, using the Python script [`logcombiner.py`](src/logcombiner.py):

		ls alignment_blocks_filtered/*/*.tre > mcc_trees.txt
		python3 logcombiner.py mcc_trees.txt mcc_trees.trees

* Open file [`mcc_trees.trees`](res/mcc_trees.trees) in DensiTree.

	<p align="center"><img src="img/densitree1.png" alt="DensiTree" width="600"></p>
	
<a name="divtimes"></a>
## Analyzing introgression with divergence times

XXX

* Sample 100 trees from the posterior distribution of each block.

		for i in alignment_blocks_filtered/*/*.trees
		do
			python3 logcombiner.py -b 10 -n 100 --remove-comments ${i} ${i%.trees}.100.trees
		done 

* Extract the node ages from each gene tree.

		mkdir node_ages
		for i in alignment_blocks_filtered/*/*.100.trees
		do
			block_id=`basename ${i%.100.trees}`
			echo -n "Reading node ages from file ${block_id}.100.trees..."
			Rscript get_mrca_table.r ${i} > node_ages/${block_id}.txt
			echo " done."
			echo "Wrote file node_ages/${block_id}.txt."
		done

* Use the Ruby script [`summarize_mrca_tables.rb`](src/summarize_mrca_tables.rb) to analyze all tables of node ages and to generate a single table of the mean pairwise node ages:

		ruby summarize_mrca_tables.rb node_ages mean_node_ages.txt

* Make a heatmap with individual sequences as units:

		R
		table <- read.table("mean_node_ages.txt")
		matrix <- as.matrix(table)
		col_palette <- colorRampPalette(c("#848e79","#6b857a","#537b7c","#3b727d","#2a667b","#255771","#204768","#1b385f","#162955","#111a4c","#0d133d","#0a0e2d"),space="rgb")(n = 25)
		col_breaks=seq(0,max(matrix),length=26)
		pdf("mean_node_ages.pdf", height=7, width=7)
		heatmap(matrix, Rowv=T, symm=T, scale="none", col=col_palette, breaks=col_breaks)
		dev.off()
		quit(save="no")

	<p align="center"><img src="img/mean_node_ages.png" alt="Mean node ages" width="600"></p>

* Prepare a file with a table assigning sequence IDs to species. Copy the following text into a new file named `species.txt`:

		IZC5_A	astbur
		IZC5_B	astbur
		AUE7_A	altfas
		AUE7_B	altfas
		JBD6_A	telvit
		JBD6_B	telvit
		JUH9_A	neobri
		JUH9_B	neobri
		LJC9_A	neocan
		LJC9_B	neocan
		KHA7_A	neochi
		KHA7_B	neochi
		IVE8_A	neocra
		IVE8_B	neocra
		JWH2_A	neogra
		JWH2_B	neogra
		JWG9_A	neohel
		JWG9_B	neohel
		JWH4_A	neomar
		JWH4_B	neomar
		JWH6_A	neooli
		JWH6_B	neooli
		ISB3_A	neopul
		ISB3_B	neopul
		ISA8_A	neosav
		ISA8_B	neosav
		KFD2_A	neowal
		KFD2_B	neowal

* Shrink the matrix with script [`shrink_matrix.rb`](src/shrink_matrix.rb):

		ruby shrink_matrix.rb mean_node_ages.txt species.txt mean_node_ages_species.txt

* Make a heatmap with species as units:

		R
		table <- read.table("mean_node_ages_species.txt")
		matrix <- as.matrix(table)
		col_palette <- colorRampPalette(c("#848e79","#6b857a","#537b7c","#3b727d","#2a667b","#255771","#204768","#1b385f","#162955","#111a4c","#0d133d","#0a0e2d"),space="rgb")(n = 25)
		col_breaks=seq(0,max(matrix),length=26)
		pdf("mean_node_ages_species.pdf", height=7, width=7)
		heatmap(matrix, Rowv=T, symm=T, scale="none", col=col_palette, breaks=col_breaks)
		dev.off()
		quit(save="no")

	<p align="center"><img src="img/mean_node_ages_species.png" alt="Mean node ages" width="600"></p>

* Calculate the maximum age reduction with the Ruby script [`check_treeness.rb`](src/check_treeness.rb).

		ruby check_treeness.rb mean_node_ages_species.txt age_reductions.txt

* Plot the age reduction:

		R
		table <- read.table("age_reductions.txt")
		matrix <- as.matrix(table)
		max(matrix)
		col_palette <- colorRampPalette(c("#848e79","#6b857a","#537b7c","#3b727d","#2a667b","#255771","#204768","#1b385f","#162955","#111a4c","#0d133d","#0a0e2d"),space="rgb")(n = 25)
		col_breaks=seq(0,max(matrix),length=26)
		pdf("age_reductions.pdf", height=7, width=7)
		heatmap(matrix, Rowv=NA, Colv=NA, symm=F, scale="none", col=col_palette, breaks=col_breaks, xlab="Taxon C (donor)", ylab="Taxon A (recipient)", margins = c(6, 6))
		dev.off()
		quit(save="no")

	<p align="center"><img src="img/age_reductions.png" alt="Age reductions" width="600"></p>
	
	<!--The upper dark row for telvit comes from comparisons with neocan. For example, in the comparison neocan (B), telvit (A), and neowal (C), the mean age for B-A is 3.005, that for B-C is 3.78, and that for A-C is 3.07.-->


* Calculate once again the maximum age reduction, after excluding "neocan" with the Ruby script [`check_treeness.rb`](src/check_treeness.rb).

		ruby check_treeness.rb mean_node_ages_species.txt age_reductions_sub1.txt neocan

* Make a new plot of the age reduction, now without "neocan":

		R
		table <- read.table("age_reductions_sub1.txt")
		matrix <- as.matrix(table)
		max(matrix)
		col_palette <- colorRampPalette(c("#848e79","#6b857a","#537b7c","#3b727d","#2a667b","#255771","#204768","#1b385f","#162955","#111a4c","#0d133d","#0a0e2d"),space="rgb")(n = 25)
		col_breaks=seq(0,max(matrix),length=26)
		pdf("age_reductions_sub1.pdf", height=7, width=7)
		heatmap(matrix, Rowv=NA, Colv=NA, symm=F, scale="none", col=col_palette, breaks=col_breaks, xlab="Taxon C (donor)", ylab="Taxon A (recipient)", margins = c(6, 6))
		dev.off()
		quit(save="no")

	<p align="center"><img src="img/age_reductions_sub1.png" alt="Age reductions" width="600"></p>

* Calculate age reductions again, excluding both "neocan" and "neopul"

		ruby check_treeness.rb mean_node_ages_species.txt age_reductions_sub2.txt neocan,neopul

* Make a new plot of the age reduction, now without "neocan" and "neopul":

		R
		table <- read.table("age_reductions_sub2.txt")
		matrix <- as.matrix(table)
		max(matrix)
		col_palette <- colorRampPalette(c("#848e79","#6b857a","#537b7c","#3b727d","#2a667b","#255771","#204768","#1b385f","#162955","#111a4c","#0d133d","#0a0e2d"),space="rgb")(n = 25)
		col_breaks=seq(0,max(matrix),length=26)
		pdf("age_reductions_sub2.pdf", height=7, width=7)
		heatmap(matrix, Rowv=NA, Colv=NA, symm=F, scale="none", col=col_palette, breaks=col_breaks, xlab="Taxon C (donor)", ylab="Taxon A (recipient)", margins = c(6, 6))
		dev.off()
		quit(save="no")

	<p align="center"><img src="img/age_reductions_sub2.png" alt="Age reductions" width="600"></p>

* Write the following text assigning sample IDs to species IDs for the three species "neobri", "neooli", and "neopul" to a new file named `species_trio1.txt`:

		JUH9_A	neobri
		JUH9_B	neobri
		JWH6_A	neooli
		JWH6_B	neooli
		ISB3_A	neopul
		ISB3_B	neopul

* Count topologies for "neobri", "neooli", "neopul":

		python3 get_topologies.py -t alignment_blocks_filtered/*/*.100.trees -s species_trio1.txt
		
* Repeat the same with another species trio.

<a name="cgenie"></a>
## Simulating recombination with c-genie

XXX

* Download c-genie recombination with c-genie:

		wget https://raw.githubusercontent.com/mmatschiner/c-genie/master/c-genie

* Make the script executable:

		chmod +x c-genie

* Have a look at the available options:

		./c-genie -h

* Run c-genie with the time calibrated SNAPP phylogeny from tutorial [Divergence-Time Estimation with SNP Data](../divergence_time_estimation_with_snp_data/README.md):

		./c-genie snapp.tre lamprologini

<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** The filtering of alignment blocks based on the proportion of missing data should have removed about 200 alignment files from directory `alignment_blocks`. In my analysis, there were 672 alignment files before filtering and 477 files remained after removing those that had more than 90% missing data.