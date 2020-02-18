#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $reftree = 'RAxML_bestTree.649KS_sequences_hmmalign_raxml_renamed.tre';

my %leaf = ();
my $in = new Bio::TreeIO(-file=>$reftree, -format=>'newick');
while(my $tree = $in->next_tree){
    foreach my $l ($tree->get_leaf_nodes){
        $leaf{$l->id} = 1;
    }
}

my %fa = ();
open my $afh, '<', '649KS_sequences_031218.fasta' or die $!;
while(<$afh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	$fa{$1} = 1;
    }
}
close $afh;

foreach my $q (keys %fa){
    print "NOT IN TREE:\t$q\n" unless(exists $leaf{$q});
}
foreach my $q (keys %leaf){
    print "NOT IN FASTA:\t$q\n" unless(exists $fa{$q});
}
