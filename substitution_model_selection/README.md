# Substitution Model Selection

A tutorial on the selection of a substitution model for phylogenetic analysis

## Summary

Before running a maximum-likelihood phylogenetic analysis, the user needs to decide which free parameters should be included in the model: should a single rate be assumed for all substitutions (as in the Jukes-Cantor model of sequence evolution; [Jukes and Cantor 1969](https://www.sciencedirect.com/science/article/pii/B9781483232119500097)) or should different rates be allowed for transitions and transversions (as in the HKY model; [Hasegawa et al. 1985](https://link.springer.com/article/10.1007/BF02101694)). Or should different rates even be used for all substitutions (as in the GTR model; [Taveré 1986](http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf))? And should the frequencies of the four nucleotides (the "state frequencies") be estimated or assumed to be all equal? The optimal number of free model parameters depends on the data available and can be chosen according to criteria such as the Akaike Information Criterion (AIC; [Akaike 1974](https://ieeexplore.ieee.org/document/1100705/)) that aim to strike a balance between improvements in model fit and the number of additional parameters required for it.

## Table of contents

* [Objective](#objective)
* [Dataset](#dataset)
* [Requirements](#requirements)

<a name="objective"></a>
## Objective

In this tutorial, I will present how to select a substitution model for phylogenetic analysis with the software [PAUP* (Swofford 2003)](http://paup.phylosolutions.com), a popular multi-utility tool for various types of phylogenetic analyses.

<a name="dataset"></a>
## Dataset

The data used in this tutorial are filtered versions of the alignments generated for 16s and rag1 sequences in tutorial [Multiple Sequence Alignment](../multiple_sequence_alignment/README.md). More information on the origin of the dataset can be found in that tutorial. As PAUP* requires alignments in Nexus format as input, use the files [`16s_filtered.nex`](data/16s_filtered.nex) and [`rag1_filtered.nex`](rag1_filtered.nex).

<a name="requirements"></a>
## Requirements

* **PAUP\*:** Installation instructions and precompiled versions of [PAUP\*](http://paup.phylosolutions.com) are available on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). First developed in the late 1980s, this software is one of the oldest programs for phylogenetic analysis, and even despite its age, its author Dave Swofford has never released a final version. Regardless, PAUP* has long been one of the most frequently used phylogenetic programs, and it is still very commonly used. Over its long lifetime, PAUP* has accumulated over [40,000 citations](https://scholar.google.ch/citations?user=H1jbCPkAAAAJ&hl=en&oi=ao) even though there is no paper for it (not even a proper manual). Until recently, PAUP* could only be purchased from Sinauer Associates for around 100 USD. Since 2015, Dave Swofford distributes updated versions of PAUP* 4.0 for free as trial versions on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). These trial versions expire after a few months, so if you would want to use PAUP* also in the future, you might have to re-download it then. This situation is likely only temporary as development is underway for PAUP* 5, which will be distributed at least in part commercially.<br>
While the descriptions in this tutorial assume that you have installed the GUI version of PAUP* for Mac OS X or Windows, it can also be followed with the command-line version of PAUP\*. If you use this command-line version, you might need to look up the equivalent commands; after starting PAUP\*, this can always be done with PAUP\*'s help text that will be displayed if you simply type "?" and hit the Enter key. The 