#!/bin/env perl

use strict;
use warnings;

my %truth = ();
my %clade2desc = ();
open my $tfh, '<', '../annotation_mc_20180727.tsv' or die $!;
while(<$tfh>){
    chomp;
    next if($_ =~ m/^#/);
    my ($id, $clade, $desc) = split(/\t/, $_);
    $clade = 'clade_not_conserved' if($clade eq 'no_clade');
    $truth{$id} = $clade;
    $clade2desc{$clade} = $desc;
}
close $tfh;

print join("\t", 'Split', 'Query', 'Clade', 'TrueClade', 'CladeMatch', 'Desc', 'TrueDesc', 'DescMatch', 'TruthInTrain')."\n";
foreach my $sf (glob("*/*tsv")){
    my %intrain = ();
    my @p = split(/\//, $sf);
    my $train = $p[0].'/'.$p[0].'_train.faa';
    open my $tfh, '<', $train or die $!;
    while(<$tfh>){
	chomp;
	if($_ =~ m/^>(.+)/){
	    $intrain{$truth{$1}} += 1;
	}
    }
    close $tfh;
    
    my $n = $sf;
    $n =~ s/.+(\d\d)\.tsv$/$1/;
    open my $sfh, '<', $sf or die $!;
    while(<$sfh>){
	chomp;
	my ($q, $clade) = split(/\t/, $_);
	$clade = 'clade_not_conserved' if($clade eq 'no_clade');
	my ($cladematch, $descmatch) = ('N', 'N');
	$cladematch = 'Y' if($clade eq $truth{$q});
	$descmatch = 'Y' if($clade2desc{$clade} eq $clade2desc{$truth{$q}});
	my $truthintrain = 0;
	$truthintrain = $intrain{$truth{$q}} if(exists $intrain{$truth{$q}});
	print join ("\t", 'split_'.$n, $q, $clade, $truth{$q}, $cladematch, $clade2desc{$clade}, $clade2desc{$truth{$q}}, $descmatch, $truthintrain)."\n";
    }
    close $sfh;
}
