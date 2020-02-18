#!/bin/env perl

use strict;
use warnings;

my $clade = './KS_rawseq_pred_training_transATPKS.clades.tsv';
my $func = './annotation_mc_20180104.tsv';

my %cladelist = ();
my %clade2func = ();

open my $ifh, '<', $clade or die $!;
while(<$ifh>){
    chomp;
    my ($name, $clade) = split(/\t/, $_);
    $clade = 'clade_not_conserved' if($clade =~ m/\|/);
    $cladelist{$clade} += 1;
}
close $ifh;

open my $ffh, '<', $func or die $!;
while(<$ffh>){
    next if ($_ =~ m/^#/);
    chomp;
    my ($name, $clade, $fa) = split("\t", $_);
    $clade2func{$clade} = $fa;
}
close $ffh;

print join("\t", 'index', 'clade_index', 'PhyloNode_index', 'function_annotation', 'color')."\n";

my $default_color = '#bdbdbd';
my %colorof = (
    'Clade_1' => '#2ca25f'
    );

my $n=1;
foreach my $k (sort keys %clade2func){
    my $fa = $clade2func{$k};
    my $color = $default_color;
    $color = $colorof{$k} if(exists $colorof{$k});
    print join("\t", $n, $k, 'NA', $fa, $color)."\n";
    $n += 1;
}
foreach my $k (sort keys %cladelist){
    unless(exists $clade2func{$k}){
	my $fa = 'no_clade';
	my $color = $default_color;
	$color = $colorof{$k} if(exists $colorof{$k});
	print join("\t", $n, $k, 'NA', $fa, $color)."\n";
	$n += 1;
    }
}
