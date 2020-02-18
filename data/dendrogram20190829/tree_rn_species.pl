#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my %p = ();
open my $afh, '<', '../dendrogram20190514/KS_precomputed_1405_hmmalign_trimmed_renamed.annotationtable.tsv' or die $!;
while(<$afh>){
    chomp;
    next if($_ =~ m/^\s/);
    my ($pathway, @clades) = split(/\t/, $_);
    my @a = ();
    foreach my $c (@clades){
	last if($c eq 'NA');
	if($c eq 'clade_not_conserved'){
	    $c = 'NC';
	}else{
	    $c =~ s/Clade_//;
	    $c = '0'.$c if($c < 10);
	}
	push @a, $c;
    }
    $p{$pathway} = join('-', @a);
}
close $afh;

my %s = ();
foreach my $f (glob("../dendrogram20190514/species*.txt")){
    open my $bfh, '<', $f or die $!;
    while(<$bfh>){
	chomp;
	my ($p, $sp) = split(/\t/, $_);
	$sp =~ s/\s+$//;
	$sp =~ s/[\\\/\(\)\[\]\{\}\^\$\*\+\?\.,:;\'\"]//g;
	#$sp =~ s/[\[\]\(\)]:,;\"\'//g;
	$s{$p} = $sp;
    }
    close $bfh;
}

my $t = new Bio::TreeIO(-file=>'./upgma.nwk', -format=>'newick');
my $o = new Bio::TreeIO(-file=>'>upgma_species.nwk', -format=>'newick');

while(my $tree = $t->next_tree){
    foreach my $leaf ($tree->get_leaf_nodes){
	if(exists $s{$leaf->id}){
	    $leaf->id('"'.$s{$leaf->id}.'"');
	    #$leaf->id($s{$leaf->id});
	}else{
	    if($leaf->id =~ m/\|/){
		if($leaf->id =~ m/^OUTGROUP/){
		    $leaf->id('"'.$leaf->id.'"');
		}else{
		    die "KEY ERROR: ".$leaf->id."\n";
		}
	    }else{
		$leaf->id('"'.$leaf->id.'"'); # reference seqs
	    }
	}
    }
    $o->write_tree($tree);
}
