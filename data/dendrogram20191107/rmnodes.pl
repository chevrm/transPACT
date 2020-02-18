#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

my %torm = ();
$torm{'Pseudomonas sp 2663|FR901464|HM047288'} = 1;
$torm{'Streptomyces himastatinicus ATCC 53653|9-methylstreptimidone|FR878059'} = 1;

my $t= new Bio::TreeIO(-file=>'upgma_species.nwk', -format=>'newick');
my $o= new Bio::TreeIO(-file=>'>upgma_species_fix.nwk', -format=>'newick');

while(my $tree = $t->next_tree){
    foreach my $leaf ($tree->get_leaf_nodes){
	if(exists $torm{$leaf->id}){
	    $tree->remove_Node($leaf);
	}else{
	    $leaf->id('"'.$leaf->id.'"');
	}
    }
    $o->write_tree($tree);
}
