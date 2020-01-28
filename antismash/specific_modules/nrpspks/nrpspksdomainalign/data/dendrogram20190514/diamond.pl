#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $protein_domain_source = 'KS_precomputed_1405_hmmalign_trimmed_renamed.fasta';

## Get pairwise PIDs for each domain
system("diamond makedb --in $protein_domain_source -d full") unless(-e 'full.dmnd');
system("diamond blastp -q $protein_domain_source -d full -e 1 -f tab -o full.dbp -p 38 -k 1000000") unless(-e 'full.dbp');
