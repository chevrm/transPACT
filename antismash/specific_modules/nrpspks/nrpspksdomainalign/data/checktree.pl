#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::SeqIO;

my $tree = new Bio::TreeIO(-file=>shift, -format=>'newick');
my $aln = new Bio::SeqIO(-file=>shift, -format=>'fasta');

my %a = ();
while(my $seq = $aln->next_seq){
    $a{$seq->id} = 1;
}

my %t = ();
while(my $t = $tree->next_tree){
    foreach my $leaf ($t->get_leaf_nodes){
	$t{$leaf->id} = 1;
    }
}

my $match = 0;
foreach my $l (keys %t){
    if(exists $a{$l}){
	$match += 1;
    }else{
	print "NO MATCH:\t$l\n";
    }
}

print "\nTREE: ".scalar(keys %t)."\nALN:  ".scalar(keys %a)."\nMATCHES: $match\n";
