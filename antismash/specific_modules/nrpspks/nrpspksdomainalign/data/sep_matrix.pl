#!/bin/env perl

use strict;
use warnings;

my $original_training_matrix = 'annotateKS_perCompound_v4_uniPK.txt';
my %intraining = ();
open my $tfh, '<', $original_training_matrix or die $!;
while(<$tfh>){
    chomp;
    unless($_ =~ m/^\s/){
	my ($cluster, @rest) = split(/\s+/, $_);
	$intraining{$cluster} = 1;
    }
}
close $tfh;

my $full_matrix = 'annotateKS_per_pred_transATPKS_mc_manualderep.txt';
my ($training_new, $novel_new) = ('annotateKS_training.txt', 'annotateKS_novel.txt');
open my $tra, '>', $training_new or die $!;
open my $nov, '>', $novel_new or die $!;
open my $ffh, '<', $full_matrix or die $!;
while(<$ffh>){
    chomp;
    if($_ =~ m/^\s/){
	print $tra "$_\n";
	print $nov "$_\n";
    }else{
	my ($cluster, @rest) = split(/\s+/, $_);
	if(exists $intraining{$cluster}){
	    print $tra "$_\n";
	}else{
	    print $nov "$_\n";
	}
    }
}
close $ffh;
close $tra;
close $nov;
