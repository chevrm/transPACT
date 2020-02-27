#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $o= new Bio::TreeIO(-file=>'upgma_trim.nwk', -format=>'newick');
my %l = ();
while(my $tree = $o->next_tree){
    foreach my $leaf ($tree->get_leaf_nodes){
	$l{$leaf->id} = 'No';
    }
}

open my $cfh, '<', 'contig_border_table.txt' or die $!;
while(<$cfh>){
    chomp;
    my ($p, $left, $right) = split(/\t/, $_);
    next if($left eq '?' || $right eq '?');
    if($left==1 || $right==1){
	my @paths = split(/_/, $p);
	foreach my $path (@paths){
	    if(exists $l{$path}){
		$l{$path} = 'Yes';
	    }else{
		print STDERR "$path\tNOT FOUND\n";
	    }
	}
    }
}
close $cfh;

foreach my $leaf (keys %l){
    print join("\t", $leaf, $l{$leaf})."\n";
}
