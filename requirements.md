# Requirements

A list of the installations required for the tutorials


## Basics

* **BASH:** The [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) shell environment is required for most of the tutorials. It is included natively in all Mac OS X and Linux distributions, but not in Windows systems. Thus, if your computer is running Windows, you will have to either use a client like [PuTTY](https://www.putty.org) to access a Linux computer remotely, or you could look into installing Linux on your computer in addition to Windows. On both Mac OS X and Linux, the standard way to use the BASH shell environment is through the application named "Terminal" that you should find in the Applications directory. As a first test if your BASH shell environment is working as it should, open the Terminal app and type `echo hello`. If the Terminal sends a greeting, the environment seems ok.

* **Ruby 2.x:** The [Ruby](https://www.ruby-lang.org/en/) programming language should also be installed by default on all Mac OS X and Linux Installations. To check the version of your Ruby installation, type `ruby --version` in a Terminal window. This should output a version number to the screen. If the version number is greater than 2.0, your installation is recent enough for the scripts used in the tutorials. If you should have an older Ruby version, please update it following the instructions on the [Ruby website](https://www.ruby-lang.org/en/).

* **Python 2.7:** The [Python](https://www.python.org) programming language should also be already installed on your machine if you're a Mac OS X or Linux system. To check if this is the case, type `python --version` in a Terminal window. If the version number given on the screen is 2.7.x, your installation of Python 2 is ready to use in the tutorials. If not, please download the latest version of Python 2 from the [Python website](https://www.python.org/downloads/).

* **Python 3.x:** Unfortunately Python 2 and Python 3 have diverged to the extent that scripts written in Python 2 can not be run with Python 3 and vice versa. Thus, Python 3 needs to be installed separately, in addition to Python 2, if it is not already installed on your machine. To check whether you already have a working installation of Python 3, type `python3 --version` in a Terminal window. If this shows a version number greater than 3.x, your Python 3 installation is ready to run the tutorials. If not, you should download the latest version of Python 3 from the [Python website](https://www.python.org/downloads/).

* **R:** If you previously installed the [R programming language](https://www.r-project.org) on your computer, make sure that you can access it through the command line by typing `R --version`. If R does not appear to be properly installed, you will find installation instructions on the [website of the R project](https://www.r-project.org).


## Libraries

* **Libraries for Python 2.7:** The [ete toolkit](http://etetoolkit.org) ([Huerta-Cepas et al. 2016](https://academic.oup.com/mbe/article/33/6/1635/2579822)) will be required for comparisons between phylogenetic trees and for topology-weighting analyses with [TWISST](https://github.com/simonhmartin/twisst) ([Martin and van Belleghem 2017](http://www.genetics.org/content/206/1/429)). Instructions for the installation of the ete toolkit on Mac OS X and Linux are provided on the [ete download webpage](http://etetoolkit.org/download/); however, the easiest way to install the ete3 toolkit might be with the pip package manager for Python, using the following command:

		python -m pip install --user ete3
		
	To ensure that the installation worked, you could execute the following command:
	
		python -c 'import ete3'
		
	If no error message is given, the ete3 library is correctly installed and ready to be used.

	A second Python library that will be required not only for analyses with TWISST but also in other tutorials is the [numpy](http://www.numpy.org) library. It can be installed in a similar way with the following command:
	
		python -m pip install --user numpy
		
	Again, the installation can be tested with this command:
	
		python -c 'import numpy'

		
* **Libraries for Python 3.x:** Libraries for Python v.3.x need to be installed separately. The required libraries are [numpy](http://www.numpy.org), [scipy](https://www.scipy.org), and [msprime](https://msprime.readthedocs.io/en/stable/index.html) [(Kelleher et al. 2016)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004842). These can be installed with pip for Python 3, using the following commands:

		python3 -m pip install --user numpy
		python3 -m pip install --user scipy
		python3 -m pip install --user msprime
		
	And the installations can be tested with these commands:
	
		python3 -c 'import numpy'
		python3 -c 'import scipy'
		python3 -c 'import msprime'


## Programs

* **MAFFT:** Installation instructions and precompiled versions of [MAFFT](https://mafft.cbrc.jp/alignment/software/) are available on [https://mafft.cbrc.jp/alignment/software/](https://mafft.cbrc.jp/alignment/software/). The installation of the program should be easy on all operating systems.

* **AliView:** To visualize sequence alignments, the software [AliView](http://www.ormbunkar.se/aliview/) ([Larsson 2014](https://academic.oup.com/bioinformatics/article/30/22/3276/2391211)) is recommended. The installation of AliView is described at [http://www.ormbunkar.se/aliview/](http://www.ormbunkar.se/aliview/) and should be possible on all operating systems.

* **BMGE:** The program [BMGE](https://research.pasteur.fr/en/software/bmge-block-mapping-and-gathering-with-entropy/) (Block Mapping and Gathering with Entropy) ([Criscuolo and Gribaldo 2010](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-210)) is highy useful to identify and remove poorly aligned regions of sequence alignments. The latest version of BMGE is provided as a Java jar file at [ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/](ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/) (choose login as guest to access the ftp server). Place this file in a convenient location on your own computer.

* **PAUP\*:** Installation instructions and precompiled versions of [PAUP\*](http://paup.phylosolutions.com) are available on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). Until recently, PAUP* could only be purchased from Sinauer Associates for around 100 USD. Since 2015, Dave Swofford distributes updated versions of PAUP* 4.0 for free as trial versions on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). These trial versions expire after a few months, so if you would want to use PAUP* also in the future, you might have to re-download it then. While the tutorial instructions will assume that you have installed the GUI version of PAUP* for Mac OS X or Windows, it can also be followed with the command-line version of PAUP\*.

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

	The tutorial instructions will assume that you compiled the parallelized PTHREADS version of RAxML (rather than the sequential or the MPI version), that you named the file simply `raxml`, and that you placed it somewhere on your computer where your system can find it (i.e. in a directory that is included in your [PATH](https://en.wikipedia.org/wiki/PATH_(variable))). One way to guarantee this on Mac OS X or Linux is to place the executable in `/usr/local/bin`, for example using (if you compiled the AVX version)
	
		mv raxmlHPC-PTHREADS-AVX /usr/local/bin/raxml
		
	To verify that the RAxML executable can be found by your system, type
	
		which raxml
		
	If this command outputs a path such as `/usr/local/bin/raxml`, the executable can be found. As another check if RAxML is working as it should, type
	
		raxml -v
		
	and you should see the version number as well as a list of contributing developers. If you do, you're ready to start the tutorial.
	
* **FigTree:** The program [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by Andrew Rambaut is a very intuitive and useful tool for the visualization and (to a limited extent) manipulation of phylogenies encoded in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format. Executables for Mac OS X, Linux, and Windows are provided on [http://tree.bio.ed.ac.uk/software/figtree/](http://tree.bio.ed.ac.uk/software/figtree/).

* **BEAST2:** The BEAST2 package, including BEAUti, BEAST2 itself, TreeAnnotator, and other tools can be downloaded from the BEAST2 website [https://www.beast2.org](https://www.beast2.org). As all these programs are written in Java, compilation is not required, and all programs should work on Mac OS X, Linux, and Windows. I recommend downloading the program versions that include the Java Runtime Environment, which may prevent conflicts with Java versions that may already be installed on your machine.<br>

* **Tracer:** Just like BEAST2, Tracer is written in Java and should work on your system without problems. The program can be downloaded for Mac OS X, Linux, or Windows from [http://tree.bio.ed.ac.uk/software/tracer/](http://tree.bio.ed.ac.uk/software/tracer/). The file with the extension `.dmg` is for Mac OS X, the one with the extension `.tgz` is for Linux, and the Windows version is the file ending in `.zip`.

* **BLAST+:** The suite of command-line tools related to BLAST (Basic Local Alignment Search Tool; [Altschul et al. 1990](https://www.sciencedirect.com/science/article/pii/S0022283605803602?via%3Dihub)) will be essential for the tutorial on ortholog detection. Installers for different platforms are available from [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (connect as anonymous guest to access these files). For installation on Mac OS X, choose the file named `ncbi-blast-2.7.1+.dmg`, for Linux, choose file `ncbi-blast-2.7.1+-1.x86_64.rpm`, and on Windows, use `ncbi-blast-2.7.1+-win64.exe`. More detailed installation instructions can be found in the [online application manual](https://www.ncbi.nlm.nih.gov/books/NBK279671/#introduction.RedHat_Linux).

* **PAML:** [PAML (Phylogenetic Analysis by Maximum Likelihood)](http://abacus.gene.ucl.ac.uk/software/paml.html) is a package of tools for phylogenetic inference developed by Ziheng Yang and colleagues ([Yang 2007](https://academic.oup.com/mbe/article/24/8/1586/1103731)). Of this package we are going to use a single tool, the program codeml, which allows rapid calculation of the [dN/dS ratio](https://en.wikipedia.org/wiki/Ka/Ks_ratio) from pairwise sequence comparisons. Downloads and installation instructions for Mac OS X, Linux, and Windows are provided on the [PAML webpage](http://abacus.gene.ucl.ac.uk/software/paml.html#download). Make sure to skip the section titled "PAML-X: A GUI for PAML" and follow the instructions in section "UNIX/Linux and Mac OSX" if you use Mac OS X or Linux, or the instructions in section "PAML for Windows 9x/NT/2000/XP/Vista/7" if you use Windows.

* **ASTRAL:** The program [ASTRAL](https://github.com/smirarab/ASTRAL) ([Zhang et al. 2017](https://link.springer.com/chapter/10.1007%2F978-3-319-67979-2_4)) allows efficient and accurate estimation of the species tree based on a set of gene trees. The latest release of ASTRAL can be obtained from [https://github.com/smirarab/ASTRAL/releases](https://github.com/smirarab/ASTRAL/releases). Click on the link for "Source code (zip)" to download the full release including the source code as well as the compiled program as a Java jar file. Within the downloaded directory you'll find a zip file with the program name and its version number. Uncompress this file and open the uncompressed directory, there you should find the ASTRAL jar file, named `astral.5.5.6.jar` or similar. Rename this file so that it is simply named `astral.jar` and place it in a convenient location on your computer.

* **aTRAM2:** The software aTRAM2 itself is easy to install as it is written in Python3. The installation is described on the [aTRAM github repository](https://github.com/juliema/aTRAM), but you can skip the part in those instructions about the virtual environment, and simply use the two commands below to download the latest version of aTRAM2 and install required Python libraries:

		git clone https://github.com/juliema/aTRAM.git
		python3 -m pip install --user -r aTRAM/requirements.txt
		
		
* **Abyss:** [Abyss](https://github.com/bcgsc/abyss) ([Jackman et al. 2017](https://genome.cshlp.org/content/27/5/768)) is one out of four interchangeable assembler programs supported and internally used by aTRAM2, and is recommended here as its installation is relatively easy on Mac OS X, Linux, and Windows. The other three supported assemblers are [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) ([Grabherr et al. 2011](https://www.nature.com/articles/nbt.1883)), [Velvet](https://www.ebi.ac.uk/%7Ezerbino/velvet/) ([Zerbino 2010](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/0471250953.bi1105s31)), and [SPAdes](http://cab.spbu.ru/software/spades/) ([Bankevich et al. 2012](https://www.ncbi.nlm.nih.gov/pubmed/22506599)); if you should have one of these already installed on your machine you could skip the installation of Abyss. The installation of of Abyss is described on the associated [github repository](https://github.com/bcgsc/abyss).

* **BWA:** The [Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net) (BWA) ([Li and Durbin 2009](https://academic.oup.com/bioinformatics/article/25/14/1754/225615)) is a requirement of Abyss. Download the latest version of BWA with the following command

		git clone https://github.com/lh3/bwa.git
		
	and then follow the instructions in `README.md` located inside the downloaded directory. Note that to compile with `make` on Mac OS X, you might need to change the first line of the file `Makefile` in the same directory from `CC = gcc` to `CC = clang` (but try the compilation without this change first).

* **bcftools:** [bcftools](http://www.htslib.org/doc/bcftools.html) ([Li 2011](https://academic.oup.com/bioinformatics/article/27/21/2987/217423)) is a fast and versatile tool for the manipulation and filtering of variant data in VCF format. Downloads and instructions for installation on Mac OS X and Linux are available at the [HTSlib download webpage](http://www.htslib.org/download/). Installation on Windows is apparently not possible. If you should fail to install bcftools, you could skip the respective tutorial part and continue with ready-made files afterwards.

* **vcftools:** Similar to bcftools, [vcftools](https://vcftools.github.io/downloads.html) ([Danecek et al. 2011](https://academic.oup.com/bioinformatics/article/27/15/2156/402296)) is a program for the manipulation of files in VCF format. vcftools is generally slower than bcftools and is increasingly being replaced by it. Nevertheless, some functions that are available in vcftools have not yet been implemented in bcftools. Download files and installation instructions can be found at the [vcftools download webpage](https://vcftools.github.io/downloads.html). Like the bcftools installation, installing vcftools is not absolutely required because you could skip the tutorial parts using it.

* **Phi Test:** [Phy Test](http://www.maths.otago.ac.nz/~dbryant/software.html) ([Bruen et al. 2006](http://www.genetics.org/content/172/4/2665)) implements a statistical test for the presence of recombination in a sequence alignment. Source code is available from [David Bryant's software webpage](http://www.maths.otago.ac.nz/~dbryant/software.html). Note that in order to compile it on a Mac, you may have to replace `CXX = gcc` with `CXX = clang` on the third line of the file `Makefile` (but try the compilation without this change first).

* **Saguaro:** [Saguaro](http://saguarogw.sourceforge.net) ([Zamani et al. 2013](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-347)) implements a hidden-Markov model for the detection of recombination breakpoints in chromosome-length alignments. Installation instructions can be found at [Saguaro's webpage](http://saguarogw.sourceforge.net), and the download is available from the associated [Sourceforge repository](https://sourceforge.net/p/saguarogw/code/HEAD/tree/) (click "Download Snapshot" to download). Note that Saguaro runs on Linux, and while in principle the installation should also be possible on a Mac, this does not seem to be easy and is not supported by the authors. The installation is not possible on Windows. If you have access to a Linux server with a Saguaro installation, but you would like to run the rest of the tutorial on your own machine, you can do so by transferring input and ouput files of Saguaro via scp between your machine and the Linux server. To check whether the Saguaro installation succeeded, just type `Saguaro` on the command line. If you should fail to install Saguaro, you can skip the tutorial part on with the Saguaro analysis and continue with ready-made Saguaro result files afterwards.

