#!/bin/env perl

use strict;
use warnings;

my %keep = ();
open my $kfh, '<', 'list2keep.txt' or die $!;
while(<$kfh>){
    chomp;
    $keep{$_} = 1;
}
close $kfh;

open my $ofh, '>', 'subset.annotationtable.tsv' or die $!;
open my $afh, '<', '../dendrogram20200320/KS_precomputed_1405_hmmalign_trimmed_renamed.annotationtable.tsv' or die $!;
while(<$afh>){
    chomp;
    if($_ =~ m/^\s/ || $_ =~ m/^OUTGROUP/){
	print $ofh "$_\n";
    }else{
	my @l = split(/\t/, $_);
	if(exists $keep{$l[0]}){
	    print $ofh join("\t", @l)."\n";
	}
    }
}
close $afh;
