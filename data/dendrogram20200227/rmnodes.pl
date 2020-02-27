#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

my %torm = ();
open my $afh, '<', '../dendrogram20190514/KS_precomputed_1405_hmmalign_trimmed_renamed.annotationtable.tsv';
while(<$afh>){
    chomp;
    next if($_ =~ m/^\s/);
    my ($pathway, @ks) = split(/\t/, $_);
    my $ks_n = 0;
    foreach my $k (@ks){
	if($k =~ m/^Clade/){
	    $ks_n += 1;
	}else{
	    last;
	}
    }
    if($ks_n <= 1){
	$torm{$pathway} = 1 unless($pathway =~ m/^OUTGROUP/);
    }
}
close $afh;

my $t= new Bio::TreeIO(-file=>'upgma.nwk', -format=>'newick');
my $o= new Bio::TreeIO(-file=>'>upgma_trim.nwk', -format=>'newick');

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
