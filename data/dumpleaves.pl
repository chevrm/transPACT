#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::SeqIO;

my $tree = new Bio::TreeIO(-file=>shift, -format=>'newick');

while(my $t = $tree->next_tree){
    foreach my $leaf ($t->get_leaf_nodes){
	print $leaf->id."\n";
    }
}
