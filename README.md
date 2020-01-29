# IMPACT
Integrative Modular trans-AT PKS Annotation and Comparison Tool


## Set up environment

`conda create --name impact --file conda_packages.txt`

## Key scripts

* `antismash/specific_modules/nrpspks/nrpspksdomainalign/substrate_from_faa.py <single_entry_fasta.faa>`
  * IMPACT prediction of trans-AT substrate from a protein fasta
* `antismash/specific_modules/nrpspks/nrpspksdomainalign/data/dendrogram20190829/generate_dendrogram_userweights.py <Jaccard_weight> <DDS_weight> <AdjacencyIndex_weight>`
  * Generate trans-AT pathway dendrogram (based on all vs all DIAMOND, not in repo)
  * Suggested weights are Jaccard = 0, DDS = 0.32, Adjacency index = 0.68

