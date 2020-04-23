# _trans_PACT
trans_-AT **P**KS **A**nnotation and **C**omparison **T**ool

<center><img src="/images/logos.png" alt="logos"
	title="transPACT is a joint collaboration between the University of Wisconsin-Madison, ETH Zurich, and Wageningen University" width="300" height="200" /></center>

transPACT is a joint collaboration between the University of Wisconsin-Madison, ETH Zurich, and Wageningen University.

### Reference:

_"Evolution of combinatorial diversity in trans-acyltransferase polyketide synthase assembly lines."_

EJN Helfrich*, R Ueoka*, MG Chevrette*, F Hemmerling, X Lu, S Leopold-Messer, AY Burch, SE Lindow, J Handelsman, J Piel†, MH Medema†

\* equal contributions

† to whom correspondance should be addressed; JP: jpiel (at) ethz.ch | MHM: marnix.medema (at) wur.nl

## Brief description

_Trans_-acyltransferase polyketide synthases (_trans_-AT PKSs) are multimodular enzymes that biosynthesize diverse pharmaceutically and ecologically important natural products. Here, we developed and applied a phylogenomic algorithm, _trans_PACT (_trans_-AT **P**KS **A**nnotation and **C**omparison **T**ool), to perform a global computational analysis of _trans_-AT PKS gene clusters, identifying hundreds of evolutionarily conserved module blocks. Network analysis of their exchange patterns reveals a widespread diversification mechanism for these enzymes. _trans_PACT implementation to assign substrate specificity to _trans_-AT PKS's ketosynthase (KS) domains can be found within this repository, as well as helper scripts used to generate the global _trans_-AT PKS network. _trans_PACT is typically run independently, but is built within the antiSMASH 4.x architecture \[[paper](https://academic.oup.com/nar/article/45/W1/W36/3778252 "Link to paper")\] \[[repo](https://bitbucket.org/antismash/antismash/src/master/ "Link to repository")\].

## Set up environment

Dependencies are listed in `conda_packages.txt`. It is highly suggested for users to create their own conda environment using this file, e.g.:

`conda create --name transPACT --file conda_packages.txt`

This creates a new environment called `transPACT` with all dependencies installed. This environment can now be accessed by:

`conda activate transPACT`

## Running _trans_PACT to assign KS substrate specificity 

* `python2 transPACT_substrate_from_faa.py <protein_fasta_of_KS_domains.faa>`
  * _trans_PACT prediction of _trans_-AT substrate from a protein fasta. An example is provided in `example/test.faa`.
  * Tab separated output (default is to STDOUT; redirect to a file to save results)
  
* `python2 ./data/dendrogram20190829/generate_dendrogram_userweights.py <Jaccard_weight> <DSS_weight> <AdjacencyIndex_weight>`
  * Generate _trans_-AT pathway dendrogram
  * Implementation of Jaccard index (JI), domain sequence similariry (DSS), and adjacency index is as described in BiG-SCAPE \[[paper](https://www.nature.com/articles/s41589-019-0400-9 "Link to paper")\]. Briefly, JI measures the percentage of shared types of domains, DSS measures sequence identity between protein domains, and AI measures the percentage of pairs of adjacent domains.
    * Suggested weights are JI = 0, DSS = 0.32, AI = 0.68, the same weights that are used in BiG-SCAPE's distance calculation for _trans_-AT PKS pathways.
  * Not provided in this repo (due to size): all vs all diamond table (filename set at line 576).
  * Output is a newick format dendrogram that can be visualized in any number of tree visualization software. We recommend iTOL, and have used that for out global analysis \[[iTOL](https://itol.embl.de/ "Link to iTOL")\] \[[our analysis](https://itol.embl.de/tree/12810415342393521582825923# "Link to our global trans-AT dendrogram")\]. Extensive documentation on annotating iTOL trees can be found [here](https://itol.embl.de/help.cgi#datasets "Link to iTOL documentation"). Our annotation files are at `data/dendrogram20200227/itol_bin.txt` for denoting whether a BGC lies on a contig edge and `data/dendrogram20200227/itol_dom.txt` for annotating the KS-domain clades of the pathway.
  
## What's actually happening when I run _trans_PACT?

The core _trans_PACT algorithm is found at `antismash/specific_modules/nrpspks/nrpspksdomainalign/substrate_from_faa.py`. It has been symbolically linked at `transPACT_substrate_from_faa.py` for user convenience. For each ketosynthase domain (input as a protein fasta), KSs are aligned to a reference alignment of a core set of 647 experimentally characterized KS domains with MUSCLE (see `align_ks_domains()`; invoked on line 533). This alignment is used to phylogenetically place the query sequence onto a reference phylogeny (placement with pplacer; see `run_pipeline_pplacer()`; invoked on line 534) and query sequences are assigned to a clade and functional classification based on monophyly (see `parse_pplacer()`).

