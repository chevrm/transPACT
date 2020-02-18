# IMPACT
**I**ntegrative **M**odular _trans_-AT **P**KS **A**nnotation and **C**omparison **T**ool

<center><img src="/images/logos.png" alt="logos"
	title="IMPACT is a joint collaboration between the University of Wisconsin-Madison, ETH Zurich, and Wageningen University" width="300" height="200" /></center>

IMPACT is a joint collaboration between the University of Wisconsin-Madison, ETH Zurich, and Wageningen University.

### Reference:

_"Evolution of combinatorial diversity in trans-acyltransferase polyketide synthase assembly lines."_

EJN Helfrich*, R Ueoka*, MG Chevrette*, F Hemmerling, X Lu, S Leopold-Messer, AY Burch, SE Lindow, J Handelsman, J Piel†, MH Medema†

\* equal contributions

† to whom correspondance should be addressed; JP: jpiel (at) ethz.ch | MHM: marnix.medima (at) wur.nl

## Brief description

_Trans_-acyltransferase polyketide synthases (_trans_-AT PKSs) are multimodular enzymes that biosynthesize diverse pharmaceutically and ecologically important natural products. Here, we developed and applied a phylogenomic algorithm, IMPACT (**I**ntegrative **M**odular _trans_-AT **P**KS **A**nnotation and **C**omparison **T**ool), to perform a global computational analysis of _trans_-AT PKS gene clusters, identifying hundreds of evolutionarily conserved module blocks. Network analysis of their exchange patterns reveals a widespread diversification mechanism for these enzymes. IMPACT implementation to assign substrate specificity to _trans_-AT PKS's ketosynthase (KS) domains can be found within this repository, as well as helper scripts used to generate the global _trans_-AT PKS network. IMPACT is typically run independently, but is built within the antiSMASH 4.x architecture \[[paper](https://academic.oup.com/nar/article/45/W1/W36/3778252 "Link to paper")\] \[[repo](https://bitbucket.org/antismash/antismash/src/master/ "Link to repository")\].

## Set up environment

Dependencies are listed in `conda_packages.txt`. It is highly suggested for users to create their own conda environment using this file, e.g.:

`conda create --name impact --file conda_packages.txt`

This creates a new environment called `impact` with all dependencies installed. This environment can now be accessed by:

`conda activate impact`

## Running IMPACT to assign KS substrate specificity 

* `antismash/specific_modules/nrpspks/nrpspksdomainalign/substrate_from_faa.py <single_entry_fasta.faa>`
  * IMPACT prediction of trans-AT substrate from a protein fasta
* `antismash/specific_modules/nrpspks/nrpspksdomainalign/data/dendrogram20190829/generate_dendrogram_userweights.py <Jaccard_weight> <DDS_weight> <AdjacencyIndex_weight>`
  * Generate trans-AT pathway dendrogram
  * Suggested weights are Jaccard = 0, DDS = 0.32, Adjacency index = 0.68
  * Not provided in this repo (due to size): all vs all diamond table (filename set at line 576).

