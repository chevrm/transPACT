#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::SeqIO;

my $tree = new Bio::TreeIO(-file=>shift, -format=>'newick');

my %a = ();
open my $ifh, '<', shift or die $!;
while(<$ifh>){
    next if ($_=~ m/^Leaf/);
    chomp;
    my ($ln, @rest) = split(/,/, $_);
    $a{$ln} = 1;
}
close $ifh;

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
	print "$l\n";
    }
}

#print "\nTREE: ".scalar(keys %t)."\nIN:  ".scalar(keys %a)."\nMATCHES: $match\n";
