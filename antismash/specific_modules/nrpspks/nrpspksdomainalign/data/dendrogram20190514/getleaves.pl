#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $reftree = shift; #'RAxML_bestTree.649KS_sequences_hmmalign_raxml.tre';

my %leaf = ();
my $in = new Bio::TreeIO(-file=>$reftree, -format=>'newick');
while(my $tree = $in->next_tree){
    foreach my $l ($tree->get_leaf_nodes){
        print $l->id."\n";
    }
}
