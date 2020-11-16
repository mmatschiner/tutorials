# Requirements

A list of the installations required for the tutorials


## Basics

* **BASH & Console:** The [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) shell environment is required for most of the tutorials. It is included natively in all Mac OS X and Linux distributions but is usually not installed on Windows systems. Thus, if your computer is running Windows, you may either have to install BASH separately (which seems to be possible since Windows 10), use a client like [PuTTY](https://www.putty.org) to access a Linux server remotely (one such server will be made available to participants of the 2020 Physalia Phylogenomics course), or you could look into installing Linux on your computer in addition to Windows. The BASH shell environment is used through a console program like the application named "Terminal" that you should find in the Applications directory if you use Mac OS X or Linux. Alternative console programs exist for these systems and could be used instead of Terminal. Equivalent program on Windows seem to be called "Windows Terminal" or "Windows Console". As a first test if your BASH shell environment is working as it should, open the console program, type `echo hello`, and hit Enter. If the program sends a greeting, the environment seems ok.

* **Ruby 2.x:** The [Ruby](https://www.ruby-lang.org/en/) programming language should also be installed by default on all Mac OS X and Linux Installations. To check the version of your Ruby installation, type `ruby --version` in a Terminal window. This should output a version number to the screen. If the version number is greater than 2.0, your installation is recent enough for the scripts used in the tutorials. If you should have an older Ruby version, please update it following the instructions on the [Ruby website](https://www.ruby-lang.org/en/).

* **Python 3.x:** The [Python](https://www.python.org) programming language should also be already installed on your machine if you're a Mac OS X or Linux system. To check if this is the case, type `python3 --version` in a Terminal window. If this shows a version number (which should be greater than 3.x), your Python 3 installation is ready to run the tutorials. If not, try whether `python --version` results in a version number greater than 3.x. If this is not the case, you should download the latest version of Python 3 from the [Python website](https://www.python.org/downloads/).

* **R:** If you previously installed the [R programming language](https://www.r-project.org) on your computer, make sure that you can access it through the command line by typing `R --version`. If R does not appear to be properly installed, you will find installation instructions on the [website of the R project](https://www.r-project.org).


## Libraries

* **Libraries for Python 3.x:** Some libraries for Python v.3.x may need to be installed separately. The required libraries are [numpy](http://www.numpy.org), [scipy](https://www.scipy.org), and [msprime](https://msprime.readthedocs.io/en/stable/index.html) [(Kelleher et al. 2016)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004842). These can be installed with pip for Python 3, using the following commands:

		python3 -m pip install --user numpy
		python3 -m pip install --user scipy
		python3 -m pip install --user msprime
		
	And the installations can be tested with these commands:
	
		python3 -c 'import numpy'
		python3 -c 'import scipy'
		python3 -c 'import msprime'


* **Libraries for R:** Two R libraries will be required for analyses of introgression; these are [ape](https://cran.r-project.org/web/packages/ape/index.html) ([Paradis et al. 2004](https://academic.oup.com/bioinformatics/article/20/2/289/204981)) and [coda](https://cran.r-project.org/web/packages/coda/index.html) ([Plummer et al. 2006](http://oro.open.ac.uk/22547/)). If you're already familiar with R, just install these packages in the usual way. If not, the easiest way to do so might be via the command line. Type `R` to open the R environment interactively. Then, run the following commands:

		install.packages("ape", repos="http://ftp.gwdg.de/pub/misc/cran/", dependencies=T)
		install.packages("coda", repos="http://ftp.gwdg.de/pub/misc/cran/", dependencies=T)
		
	To ensure that both packages have been successfully installed, type these commands:
	
		library(ape)
		library(coda)

	If both commands result in no error messages, the two packages are ready to be used. Then, quit the R environment with `quit(save="no")`.


## Programs

* **MAFFT:** Installation instructions and precompiled versions of [MAFFT](https://mafft.cbrc.jp/alignment/software/) are available on [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/). The installation of the program should be easy on all operating systems.

* **AliView:** To visualize sequence alignments, the software [AliView](http://www.ormbunkar.se/aliview/) ([Larsson 2014](https://academic.oup.com/bioinformatics/article/30/22/3276/2391211)) is recommended. The installation of AliView is described at [http://www.ormbunkar.se/aliview/](http://www.ormbunkar.se/aliview/) and should be possible on all operating systems.

* **BMGE:** The program [BMGE](https://research.pasteur.fr/en/software/bmge-block-mapping-and-gathering-with-entropy/) (Block Mapping and Gathering with Entropy) ([Criscuolo and Gribaldo 2010](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-210)) is highy useful to identify and remove poorly aligned regions of sequence alignments. The latest version of BMGE is provided as a Java jar file at [ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/](ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/) (choose login as guest to access the ftp server). Place this file in a convenient location on your own computer.

* **PAUP\*:** Installation instructions and precompiled versions of [PAUP\*](http://paup.phylosolutions.com) are available on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). Until recently, PAUP* could only be purchased from Sinauer Associates for around 100 USD. Since 2015, Dave Swofford distributes updated versions of PAUP* 4.0 for free as trial versions on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). These trial versions expire after a few months, so if you would want to use PAUP* also in the future, you might have to re-download it then. While the tutorial instructions will assume that you have installed the GUI version of PAUP* for Mac OS X or Windows, it can also be followed with the command-line version of PAUP\*.

* **IQ-TREE:** The IQ-TREE version 2.0.6 or greater will be required as earlier versions did not allow the calculations of gene and site concordance factors. The latest versions of IQ-TREE for Mac OS X, Linux, or Windows can be found on [https://github.com/Cibiv/IQ-TREE/releases](https://github.com/Cibiv/IQ-TREE/releases). To install IQ-TREE on any of these systems, download the version for your operating system, and decompress this file on your machine if necessary. In the decompressed directory, you'll find a subdirectory named `bin` and inside of this subdirectory should be a file named `iqtree` or `iqtree.exe`.
	
* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by Andrew Rambaut is a very intuitive and useful tool for the visualization and (to a limited extent) manipulation of phylogenies encoded in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format. Executables for Mac OS X, Linux, and Windows are provided on [https://github.com/rambaut/figtree/releases](https://github.com/rambaut/figtree/releases).

* **BEAST2:** The BEAST2 package, including BEAUti, BEAST2 itself, TreeAnnotator, and other tools can be downloaded from the BEAST2 website [https://www.beast2.org](https://www.beast2.org). As all these programs are written in Java, compilation is not required, and all programs should work on Mac OS X, Linux, and Windows. I recommend downloading, where possible, the program versions that include the Java Runtime Environment, which may prevent conflicts with Java versions that may already be installed on your machine.<br>

* **Tracer:** Just like BEAST2, Tracer is written in Java and should work on your system without problems. The program can be downloaded for Mac OS X, Linux, or Windows from [https://github.com/beast-dev/tracer/releases](https://github.com/beast-dev/tracer/releases). The file with the extension `.dmg` is for Mac OS X, the one with the extension `.tgz` is for Linux, and the Windows version is the file ending in `.zip`.

* **BLAST+:** The suite of command-line tools related to BLAST (Basic Local Alignment Search Tool; [Altschul et al. 1990](https://www.sciencedirect.com/science/article/pii/S0022283605803602?via%3Dihub)) will be essential for the tutorial on ortholog detection. Installers for different platforms are available from [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (connect as anonymous guest to access these files). For installation on Mac OS X, choose the file named `ncbi-blast-2.11.0+.dmg`, for Linux, choose file `ncbi-blast-2.11.0+-1.x86_64.rpm`, and on Windows, use `ncbi-blast-2.11.0+-win64.exe`. More detailed installation instructions can be found in the [online application manual](https://www.ncbi.nlm.nih.gov/books/NBK279671/#introduction.RedHat_Linux).

* **PAML:** [PAML (Phylogenetic Analysis by Maximum Likelihood)](http://abacus.gene.ucl.ac.uk/software/paml.html) is a package of tools for phylogenetic inference developed by Ziheng Yang and colleagues ([Yang 2007](https://academic.oup.com/mbe/article/24/8/1586/1103731)). Of this package we are going to use a single tool, the program codeml, which allows rapid calculation of the [dN/dS ratio](https://en.wikipedia.org/wiki/Ka/Ks_ratio) from pairwise sequence comparisons. Downloads and installation instructions for Mac OS X, Linux, and Windows are provided on the [PAML webpage](http://abacus.gene.ucl.ac.uk/software/paml.html#download). Make sure to skip the section titled "PAML-X: A GUI for PAML" and follow the instructions in section "UNIX/Linux and Mac OSX" if you use Mac OS X or Linux, or the instructions in section "PAML for Windows 9x/NT/2000/XP/Vista/7" if you use Windows.

* **ASTRAL:** The program [ASTRAL](https://github.com/smirarab/ASTRAL) ([Zhang et al. 2017](https://link.springer.com/chapter/10.1007%2F978-3-319-67979-2_4)) allows efficient and accurate estimation of the species tree based on a set of gene trees. The latest release of ASTRAL can be obtained from [https://github.com/smirarab/ASTRAL/releases](https://github.com/smirarab/ASTRAL/releases). Click on the link for "Source code (zip)" to download the full release including the source code as well as the compiled program as a Java jar file. Within the downloaded directory you'll find a zip file with the program name and its version number. Uncompress this file and open the uncompressed directory, there you should find the ASTRAL jar file, named `astral.5.7.1.jar` or similar. Rename this file so that it is simply named `astral.jar` and place it in a convenient location on your computer.

* **aTRAM2:** The software aTRAM2 itself is easy to install as it is written in Python3. The installation is described on the [aTRAM github repository](https://github.com/juliema/aTRAM), but you can skip the part in those instructions about the virtual environment, and simply use the two commands below to download the latest version of aTRAM2 and install required Python libraries:

		git clone https://github.com/juliema/aTRAM.git
		python3 -m pip install --user -r aTRAM/requirements.txt
		
		
* **Abyss:** [Abyss](https://github.com/bcgsc/abyss) ([Jackman et al. 2017](https://genome.cshlp.org/content/27/5/768)) is one out of four interchangeable assembler programs supported and internally used by aTRAM2, and is recommended here as its installation is relatively easy on Mac OS X, Linux, and Windows. The other three supported assemblers are [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) ([Grabherr et al. 2011](https://www.nature.com/articles/nbt.1883)), [Velvet](https://www.ebi.ac.uk/%7Ezerbino/velvet/) ([Zerbino 2010](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/0471250953.bi1105s31)), and [SPAdes](http://cab.spbu.ru/software/spades/) ([Bankevich et al. 2012](https://www.ncbi.nlm.nih.gov/pubmed/22506599)); if you should have one of these already installed on your machine you could skip the installation of Abyss. The installation of of Abyss is described on the associated [github repository](https://github.com/bcgsc/abyss).

* **BWA:** The [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net) (BWA) ([Li and Durbin 2009](https://academic.oup.com/bioinformatics/article/25/14/1754/225615)) is a requirement of Abyss. Download the latest version of BWA with the following command

		git clone https://github.com/lh3/bwa.git
		
	and then follow the instructions in `README.md` located inside the downloaded directory. Note that to compile with `make` on Mac OS X, you might need to change the first line of the file `Makefile` in the same directory from `CC = gcc` to `CC = clang` (but try the compilation without this change first).

* **bcftools:** [bcftools](http://www.htslib.org/doc/bcftools.html) ([Li 2011](https://academic.oup.com/bioinformatics/article/27/21/2987/217423)) is a fast and versatile tool for the manipulation and filtering of variant data in VCF format. Downloads and instructions for installation on Mac OS X and Linux are available at the [HTSlib download webpage](http://www.htslib.org/download/). Installation on Windows is apparently not possible. If you should fail to install bcftools, you could skip the respective tutorial part and continue with ready-made files afterwards.

* **Dsuite:** The [Dsuite](https://github.com/millanek/Dsuite) program allows the fast calculation of the *D*-statistic from SNP data in VCF format. The program is particularly useful because it automatically calculates the *D*-statistic either for all possible species quartets or for subsets of quartets that are compatible with a user-provided species tree. Instructions for download and installation on Mac OS X and Linux are provided on [https://github.com/millanek/Dsuite](https://github.com/millanek/Dsuite). Installation on Windows is not supported, but Windows users could use ready-made output files (which will be provided) to learn how to plot and analyze the Dsuite output.
 
* **Relate:** The software [Relate](https://myersgroup.github.io/relate/index.html) estimates genome-wide sets of genealogies and is extremely fast in doing so. Downloads for Mac OS X and Linux are provided on [https://myersgroup.github.io/relate/index.html](https://myersgroup.github.io/relate/index.html), installation on Windows is apparently not supported. After downloading and decompressing the file containing the software, a collection of compiled programs (named `Relate`, `RelateMutationRate` etc.) can be found in the `bin` subdirectory inside of the decompressed directory.