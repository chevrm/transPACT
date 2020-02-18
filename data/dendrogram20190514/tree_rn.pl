#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $t = new Bio::TreeIO(-file=>'./RAxML_bestTree.649KS_sequences_hmmalign_raxml.tre', -format=>'newick');
my $o = new Bio::TreeIO(-file=>'>RAxML_bestTree.649KS_sequences_hmmalign_raxml_renamed.tre', -format=>'newick');

while(my $tree = $t->next_tree){
    foreach my $leaf ($tree->get_leaf_nodes){
	my $rn = $leaf->id;
	$rn =~ s/\.//g;
	$leaf->id($rn);
    }
    $o->write_tree($tree);
}
