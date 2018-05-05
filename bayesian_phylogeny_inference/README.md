# Bayesian Phylogenetic Inference

A tutorial on Bayesian inference of time-calibrated phylogenies

## Summary

XXX

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [Bayesian phylogenetic inference with BEAST2](#beast2)


<a name="outline"></a>
## Outline

In this tutorial, I will demonstrate how time-calibrated phylogenies can be inferred with programs of the Bayesian software package [BEAST2](https://www.beast2.org) ([Bouckaert et al. 2014](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003537)). The settings for the analysis will be specified with the program BEAUti, the Bayesian analysis itself is going to be conducted with BEAST2, and a summary tree will be generated with the program TreeAnnotator. These three programs are part of the BEAST2 package. In addition, I will present the use of the program [Tracer](http://tree.bio.ed.ac.uk/software/tracer/) ([Rambaut et al. 2018](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy032/4989127)) to assess stationarity of the Bayesian analysis.


<a name="dataset"></a>
## Dataset

The data used in this tutorial are the filtered versions of the alignments generated for 16s and rag1 sequences in tutorial [Multiple Sequence Alignment](../multiple_sequence_alignment/README.md). These alignments contain sequence data for 41 teleost fish species and are 486 and 1,368 bp long, respectively. More information on the origin of the dataset can be found in the [Multiple Sequence Alignment](../multiple_sequence_alignment/README.md) tutorial. The software BEAST2 requires input in Nexus format, therefore we will use the files [`16s_filtered.nex`](data/16s_filtered.nex) and [`rag1_filtered.nex`](data/rag1_filtered.nex).

<a name="requirements"></a>
## Requirements

* **BEAST2:** The BEAST2 package, including BEAUti, BEAST2 itself, TreeAnnotator, and other tools can be downloaded from the BEAST2 website [https://www.beast2.org](https://www.beast2.org). As all these programs are written in Java, compilation is not required, and all programs should work on Mac OS X, Linux, and Windows. I recommend downloading the program versions that include the Java Runtime Environment, which may prevent conflicts with Java versions that may already be installed on your machine.<br>

* **Tracer:** Just like BEAST2, Tracer is written in Java and should work on your system without problems. The program can be downloaded for Mac OS X, Linux, or Windows from [http://tree.bio.ed.ac.uk/software/tracer](http://tree.bio.ed.ac.uk/software/tracer).

<a name="beast2"></a>
## Bayesian phylogenetic inference with BEAST2

In this part of the tutorial, we will run a basic Bayesian phylogenetic analysis for the 16s and rag1 alignments, using the programs of the BEAST2 software package. In this analysis, we are going to assume that the 16s and rag1 genes share the same evolutionary history and that this history is also identical to the evolutionary history of the 41 teleost fish species from which the sequences were obtained. In other words, we are going to assume that the two "gene trees" are identical and that they also are identical to the "species tree".

* If you're not familiar yet with Bayesian analyses and Markov-Chain Monte Carlo methods in general, you might be overwhelmed at first by the complexities of this type of analyses. Thus, it might be worth noting the many resources made available by BEAST2 authors that provide a wealth of information, that, even if you don't need them right now, could prove to be useful at a later stage. Thus, you might want to take a moment to explore the [BEAST2 website](https://www.beast2.org) and quickly browse through the [glossary of terms related to BEAST2 analyses](https://www.beast2.org/glossary/index.html). Note that the BEAST2 website also provides a [wide range of tutorials and manuals](https://www.beast2.org/tutorials/index.html). In addition, you can find many further tutorials on the [Taming the BEAST](https://taming-the-beast.org) website, where you will also find information about the excellent [Taming-the-BEAST workshops](https://taming-the-beast.org/workshops/). Finally, if you have further questions regarding BEAST2, you could have a look if somebody else already asked those questions on the very active [user forum](https://groups.google.com/forum/#!forum/beast-users), or you could ask these questions there yourself.

* Open the program BEAUti from the BEAST2 package, and import both the filtered 16s and rag1 alignments. To do so, click "Import Alignment" from the "File" menu and select the two files [`16s_filtered.nex`](data/16s_filtered.nex) and [`rag1_filtered.nex`](data/rag1_filtered.nex). The BEAUti window should then look as shown in the screenshot below.<p align="center"><img src="img/beauti1.png" alt="BEAUti" width="700"></p>