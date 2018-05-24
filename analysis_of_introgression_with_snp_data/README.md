# Analysis of Introgression with SNP Data

A tutorial on the analysis of hybridization and introgression with SNP data

## Summary

XXX

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Assessing introgression with D-statistics](#dstatistics)
* [Allele phasing](#phasing)
* [Topology weighting with TWISST](#twisst)
* [Using alleles fixed in parental species](#fixedalleles)
* [Improved allele phasing for recent hybrids](#improvedphasing)


<a name="outline"></a>
## Outline

In this tutorial I am going XXX


<a name="dataset"></a>
## Dataset

XXX

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

* **BEAGLE:** [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) ([Browning and Browning 2016](https://www.cell.com/ajhg/fulltext/S0002-9297(15)00491-7))

* **Python library ete3:** XXX
		

<a name="dstatistics"></a>
## Assessing introgression with D-statistics
		
* Download Simon Martin's "genomics_general" github repository.

		git clone https://github.com/simonhmartin/genomics_general.git
		
* Add it to the Python path so that its functions can be found by other Python scripts:

		export PYTHONPATH=genomics_general
		
* Convert the VCF file into genotype file format.

		python genomics_general/VCF_processing/parseVCF.py -i NC_031969.f5.sub1.vcf.gz | gzip > NC_031969.f5.sub1.geno.gz

* Write the following block to a new file named `samples.txt`:

		IZA1	astbur
		IZC5	astbur
		AUE7	altfas
		AXD5	altfas
		JBD5	telvit
		JBD6	telvit
		JUH9	neobri
		JUI1	neobri
		LJC9	neocan
		LJD1	neocan
		KHA7	neochi
		KHA9	neochi
		IVE8	neocra
		IVF1	neocra
		JWH1	neogra
		JWH2	neogra
		JWG8	neohel
		JWG9	neohel
		JWH3	neomar
		JWH4	neomar
		JWH5	neooli
		JWH6	neooli
		ISA6	neopul
		ISB3	neopul
		ISA8	neosav
		IYA4	neosav
		KFD2	neowal
		KFD4	neowal


* Define the population IDs.

		spc1=altfas
		spc2=neocan
		spc3=telvit
		spc4=astbur

* Calculate allele frequencies.

		python genomics_general/freq.py -g NC_031969.f5.sub1.geno.gz -p ${spc1} -p ${spc2} -p ${spc3} -p ${spc4} --popsFile samples.txt --target derived | grep -v nan | gzip > NC_031969.f5.sub1.pops1.tsv.gz
		

* Write a file with the chromosome lengths. Name it `chr_lengths.txt`:

		NC_031965	38372991
		NC_031966	35256741
		NC_031967	14041792
		NC_031968	54508961
		NC_031969	38038224
		NC_031970	34628617
		NC_031971	44571662
		NC_031972	62059223
		NC_031973	30802437
		NC_031974	27519051
		NC_031975	32426571
		NC_031976	36466354
		NC_031977	41232431
		NC_031978	32337344
		NC_031979	39264731
		NC_031980	36154882
		NC_031981	40919683
		NC_031982	37007722
		NC_031983	31245232
		NC_031984	36767035
		NC_031985	37011614
		NC_031986	44097196
		NC_031987	43860769
		UNPLACED	141274046



* Calculate D-statistics over the entire genome.

		Rscript calculate_abba_baba.r NC_031969.f5.sub1.pops1.tsv.gz pops1.abba_baba.txt ${spc1} ${spc2} ${spc3} ${spc4} chr_lengths.txt

* Calculate D-statistics in sliding windows.

		python genomics_general/ABBABABAwindows.py -g NC_031969.f5.sub1.geno.gz -o pops1.abba_baba.windows.txt -P1 ${spc1} -P2 ${spc2} -P3 ${spc3} -O ${spc4} --popsFile samples.txt -f phased -w 500000 --windType coordinate --minSites 500

* Plot D-statistics in sliding windows.

		Rscript plot_abbababa_windows.r pops1.abba_baba.windows.txt chr_lengths.txt pops1.abba_baba.windows.pdf
		
* Repeat for different sets of species, use "pops2", "pops3" etc. instead of "pops1" in the output file names.

**Question 1:** XXX? [(see answer)](#q1)

<a name="phasing"></a>
## Allele phasing

XXX

* Download the Java jar file for BEAGLE:

		wget https://faculty.washington.edu/browning/beagle/beagle.16May18.771.jar

* Run imputation with BEAGLE:

		java -jar -Xmx4G beagle.16May18.771.jar nthreads=1 chrom=NC_031969 ne=100000 gt="NC_031969.f5.sub1.vcf.gz" out="NC_031969.f5.sub1.phased"
		
	This analysis should take around 7 minutes.

* Mask the imputed genotypes.

		gunzip -c NC_031969.f5.sub1.vcf.gz | grep -v "#" > original.vcf
		gunzip -c NC_031969.f5.sub1.vcf.gz | grep "#" > header.vcf
		gunzip -c NC_031969.f5.sub1.phased.vcf.gz | grep -v "#" > phased.vcf
		ruby mask_imputed_gts.rb original.vcf phased.vcf masked.vcf
		cat header.vcf masked.vcf | gzip > NC_031969.f5.sub1.phased.masked.vcf.gz
		
* Clean up.
		
		rm original.vcf header.vcf phased.vcf masked.vcf

<a name="twisst"></a>
## Topology weighting with TWISST

* Get the github repository:

		git clone https://github.com/simonhmartin/twisst.git

* Convert the phased and masked file to geno format.

		python genomics_general/VCF_processing/parseVCF.py -i NC_031969.f5.sub1.phased.masked.vcf.gz | gzip > NC_031969.f5.sub1.phased.masked.geno.gz
		python genomics_general/filterGenotypes.py -s IZA1,IZC5,AUE7,AXD5,JBD5,JBD6,LJC9,LJD1 -i NC_031969.f5.sub1.phased.masked.geno.gz | grep -v "N|N" > pops1.geno

* Run RAxML in sliding windows:

		python genomics_general/phylo/raxml_sliding_windows.py -T 1 -g pops1.geno --outgroup IZA1,IZC5 --prefix pops1.raxml -w 100 --windType sites -f phased
		
* Write groups file, name it `pops1.samples.txt`:

		IZA1_A	astbur
		IZA1_B	astbur
		IZC5_A	astbur
		IZC5_B	astbur
		AUE7_A	altfas
		AUE7_B	altfas
		AXD5_A	altfas
		AXD5_B	altfas
		JBD5_A	telvit
		JBD5_B	telvit
		JBD6_A	telvit
		JBD6_B	telvit
		LJC9_A	neocan
		LJC9_B	neocan
		LJD1_A	neocan
		LJD1_B	neocan


* Run TWISST:

		python twisst/twisst.py -t pops1.raxml.trees.gz -w pops1.weights.csv.gz -g ${spc1} -g ${spc2} -g ${spc3} -g ${spc4} --groupsFile pops1.samples.txt --method complete

* Plot TWISST results:

		Rscript plot_twisst_per_lg.r pops1.weights.csv.gz pops1.raxml.data.tsv 38038224 pops1.samples.txt pops1.smooth.pdf pops1.rect.pdf

<a name="fixedalleles"></a>
## Using alleles fixed in parental species

* Uncompress the VCF file.

		gunzip -c NC_031969.f5.sub1.phased.masked.vcf.gz > NC_031969.f5.sub1.phased.masked.vcf

* Run the Ruby script [`get_fixed_site_gts.rb`](src/get_fixed_site_gts.rb) to determine the alleles at sites that are fixed in the two parents.

		ruby get_fixed_site_gts.rb NC_031969.f5.sub1.phased.masked.vcf pops1.fixed.txt AUE7,AXD5 JBD5,JBD6 LJC9,LJD1 1.0

* Run the Ruby script [`plot_fixed_site_gts.rb`](src/plot_fixed_site_gts.rb) to plot the genotypes of sites that are fixed in the two parents.

		ruby plot_fixed_site_gts.rb pops1.fixed.txt pops1.fixed.svg 1.0 1000


<a name="improvedphasing"></a>
## Improved allele phasing for recent hybrids

* Extract separately the header and the body of the VCF file.

		cat NC_031969.f5.sub1.phased.masked.vcf | grep "#" > header.vcf
		cat NC_031969.f5.sub1.phased.masked.vcf | grep -v "#" > main.vcf
		
* Run Ruby script [`fix_hybrid_phasing.rb`](src/fix_hybrid_phasing.rb).

		ruby fix_hybrid_phasing.rb header.vcf main.vcf main.mod.vcf samples.txt altfas,neocan,telvit

* Combine again the header and the main part of the modified VCF.

		cat header.vcf main.mod.vcf > NC_031969.f5.sub1.phased.masked.mod.vcf

* Clean up.

		rm header.vcf main.vcf main.mod.vcf

* Generate once again a plot of the genotypes at fixed sites.

		ruby get_fixed_site_gts.rb NC_031969.f5.sub1.phased.masked.mod.vcf pops1.mod.fixed.txt AUE7,AXD5 JBD5,JBD6 LJC9,LJD1 1.0
		ruby plot_fixed_site_gts.rb pops1.mod.fixed.txt pops1.mod.fixed.svg 1.0 1000


<br><hr>

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** XXX