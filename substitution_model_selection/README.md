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

* **PAUP\*:** Installation instructions and precompiled versions of [PAUP\*](http://paup.phylosolutions.com) are available on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). First developed in the late 1980s, this software is one of the oldest programs for phylogenetic analysis, and even despite its age, its author Dave Swofford has never released a final version. Regardless, PAUP* has long been one of the most frequently used phylogenetic programs, and it is still very commonly used. Over its long lifetime, PAUP* has accumulated over [40,000 citations](https://scholar.google.ch/citations?user=H1jbCPkAAAAJ&hl=en&oi=ao) even though there is no paper for it (not even a proper manual). Until recently, PAUP* could only be purchased from Sinauer Associates for around 100 USD. Since 2015, Dave Swofford distributes updated versions of PAUP* 4.0 for free as trial versions on [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/). These trial versions expire after a few months, so if you would want to use PAUP* also in the future, you might have to re-download it then. This situation is likely only temporary as development is underway for PAUP* 5, which will be distributed at least in part commercially (fun fact: While the acronym PAUP* used to stand for "Phylogenetic Analysis using Parsimony", this has recently changed to the clever recursive acronym "Phylogenetic Analysis Using PAUP").<br>
While the descriptions in this tutorial assume that you have installed the GUI version of PAUP* for Mac OS X or Windows, it can also be followed with the command-line version of PAUP\*. If you use this command-line version, you might need to look up the equivalent commands; after starting PAUP\*, this can always be done with PAUP\*'s help screen that will be displayed if you simply type "?" and hit the Enter key. The screenshot below shows the help screen of the command-line versions of PAUP\*.
<p align="center"><img src="img/paup1.png" alt="Command-line version of PAUP*" width="600"></p>

<a name="paup"></a>
## Model selection and basic phylogeny inference in PAUP*

Comparisons of substitution models based on their fit to sequence data has been implemented in several tools, and has most often been performed using the program [jModelTest](https://github.com/ddarriba/jmodeltest2) ([Darriba et al. 2012](https://www.nature.com/articles/nmeth.2109)). But since automatic selection of substitution models has recently also been implemented in the PAUP\* and the installation of PAUP\* will anyway be required for other tutorials in this repository, I here present model selection with PAUP\* instead of jModelTest. In practice, the model selection is very similar between the two programs.

* Click "Open..." in the "File" menu of PAUP\*. Make sure that at the bottom of the opening window, "Execute" is selected as the initial mode as shown in the next screenshot. Select the file with the aligned 16s sequences in Nexus format ([`16s_filtered.nex`](data/16s_filtered.nex)) and click "Open". PAUP* will give a short report of its interpretation of the file, including the number of species (taxa) and characters found in the alignment.
<p align="center"><img src="img/paup2.png" alt="PAUP*" width="600"></p>

* The option for "Automatic Model Selection" can be found in PAUP\*'s "Analysis" menu. However, when you click on it, you'll see that in order to run this model selection, a phylogeny is required. While this might appear as if it would possibly leads to circular inference (selecting a substitution model is required for maximum-likelihood phylogenetic analysis but also depends on a phylogeny), this is not an issue in practice because the outcome of the model selection does not depend strongly on having the correct phylogeny; thus, any reasonable phylogeny will lead to a similar result of the model selection. Thus, the best solution is to run a quick phylogenetic analysis with the Neighbor-Joining algorithm ([Saitou and Nei 1987](https://academic.oup.com/mbe/article/4/4/406/1029664)), which is conveniently also implemented in PAUP\*.
<p align="center"><img src="img/paup3.png" alt="PAUP*" width="600"></p>

* To choose from the available settings for Neighbor-Joining phylogenetic analysis, click on "Neighbor Joining/UPGMA..." in PAUP\*'s "Analysis" menu.
<p align="center"><img src="img/paup4.png" alt="PAUP*" width="600"></p>

* In the newly opened pop-up window, keep all default options and click "OK".

* Click once more on "Automated Model Selection..." in the "Analysis" menu. The tree generated with Neighbor-Joining will already be selected for use in the model selection, and the pop-up window will now give you several options for this model selection. The available criteria for model selection are called "AIC", "AICc", "BIC", and "DT". These are similar to likelihood-ratio tests but have the advantage that they can be used to compare models that are not "nested" (two models are nested if one of them has all the parameters of the other models plus additional parameters). "AIC" stands for the "Akaike information criterion", "AICc" is the "Akaike information criterion corrected for small sample sizes", "BIC" is the "Bayesian information criterion", and "DT" is a "decision-theoretic" criterion. The most commonly used of these is the Akaike Information Criterion ([Akaike 1974](https://ieeexplore.ieee.org/document/1100705/)). The AIC of each model is calculated independently as:<br><br>
AIC = 2 *k* &minus;2 log(*L*),<br><br>
where *k* is the number of free parameters in a model and *L* is the likelihood after all free parameters have been optimized (i.e. the maximum likelihood). Usually, a model is considered preferrable over another model if its AIC score is at least 4 points better (= smaller) than the AIC score of the other model.<br>
Set the tick next to "AIC" but remove the ticks next to "AICc", "BIC", and "DT". Also select "AIC" to the right of "Apply settings for model chosen by:". As the "Model set", choose the number "3". This means that models with equal substitution rates (e.g. the Jukes-Cantor model), with separate substitution rates for transitions and transversion (e.g. the HKY model), and models with six independent substitution rates (the GTR model) will be tested. Keep the ticks next to "equal rates" and "gamma" (allowing gamma-distributed among-site rate variation), but remove the ticks for "invar. sites" and "both". I recommend doing so because the parameters for the proportion of invariable sites ("+I") and for among-site rate-variation ("+G") are confounded as applying a particularly low rate to a set of sites has nearly the same effect as considering these sites entirely invariable. Keep the tick next to "Show output for each model", and also set the tick next to "Show parameter estimates for each model". Make sure that the settings panel looks as shown in the below screenshot, then click "OK".
<p align="center"><img src="img/paup5.png" alt="PAUP*" width="600"></p>

* PAUP\* will report the output of the model selection in three tables. In the first (under "Evaluating models for tree 1"), you'll see a list of twelve models that were compared, as shown below ("JC" stands for the Jukes-Cantor model). **Question 1:** Which model has the highest log likelihood? [(see answer)](#q1)
<p align="center"><img src="img/paup6.png" alt="PAUP*" width="600"></p>

* In columns 4 and 5 of the same table, you'll see *k*, the number of free parameters in the model. Column 4 lists the number of additional free parameters compared to the simplest model, column 5 lists the total number of free parameters. **Question 2:** Even for the simplest model, the Jukes-Cantor model, a total of 79 free parameters are listed, do you know what these could be? [(see answer)](#q2)

* The second table lists the parameter estimates for each model. After the number and the name of each model, there are nine columns with numbers. **Question 3:** Do you know what these numbers mean? [(see answer)](#q3)

* Finally, the third table lists once again the models, but this time ranked by their AIC score. **Question 4:** Which model fits best to the data? [(see answer)](#q4) **Question 5:** What is the difference in AIC score to the second-best model? [(see answer)](#q5)
<p align="center"><img src="img/paup7.png" alt="PAUP*" width="600"></p>

* Repeat the comparison of substitution models with the alignment of rag1 sequences ([rag1_filtered.nex](data/rag1_filtered.nex)). **Question 6:** Which model is the best-fitting model for this alignment, according to AIC? [(see answer)](#q6)

<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>

## Answers

<a name="q1"></a>

* **Question 1:** The GTR model with gamma-distributed rate variation has the highest log likelihood (4610.77).

<a name="q2"></a>

* **Question 2:** These are the branch lengths (scaled by the substitution rate) of the phylogeny, which are optimized to calculate the maximum likelihood for each model. Each unrooted phylogeny of *N* taxa has 2 *N* &minus; 3 branches, therefore the unrooted phylogeny for the 41 species in the 16S alignment has
2 &times; 41 &minus; 3 = 79 branches.

<a name="q3"></a>

* **Question 3:** These first four of these columns represent the stationary frequencies assumed or estimated for the four nucleotides A, C, G, and T. This is followed by six columns for the assumed or estimated rates for the substitutions A &rarr; C, A &rarr; G, A &rarr; T, C &rarr; G, C &rarr; T, and G &rarr; T. In all models these rates are assumed to be time-reversible, meaning that e.g. the substitution A &rarr; C has the same rate as C &rarr; A. All rates are measured relative to the rate of G &rarr; T substitutions, therefore the last of the six columns shows "1" for all models. The second-last column contains the values for the alpha parameter of the gamma distribution for among-site rate variation only for those models that included this parameter ("+G"). The very last column is empty because a parameter for the proportion of invariable sites was not included in any of the models.

<a name="q4"></a>

* **Question 4:** According to the AIC, the GTR model with gamma-distributed rate variation has the best fit to the data, with an AIC score of 9,397.54.

<a name="q5"></a>

* **Question 5:** The model with the second-best AIC score is the SYM model with gamma-distributed rate variation that allows six different substitution rates but assumes that the stationary frequences of all four nucleotides is equal. The difference in AIC score compared to the best-fitting model is 13.3.

<a name="q6"></a>

* **Question 6:** The GTR model with gamma-distributed rate variation has again the lowest AIC score, indicating the best fit to the sequence data. Given that even with the modest sizes of the alignments used here the complex GTR model with rate variation is selected over other models, this model is highly likely to also fit best to alignments with more taxa and longer sequences that are commonly used in phylogenomic studies.