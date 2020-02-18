#!/bin/env perl

use strict;
use warnings;


my @ko = ();
for my $n (1..38){
    push @ko, 'KS'.$n;
}
for my $n (1..4){
    push @ko, 'KS'.$n.'b';
}


my ($ann, $fasta) = (shift, shift);
my %inann = ();
open my $afh, '<', $ann or die $!;
while(<$afh>){
    chomp;
    next if($_ =~ m/^\s/);
    my ($gk, $kscount, @ks) = split(/\t/, $_);
    my $n = 0;
    foreach my $clade (@ks){
	unless($clade eq 'NA'){
	    my $ksname = $gk.'|'.$ko[$n];
	    $inann{$ksname} = $clade;
	}
	$n += 1;
    }
}
close $afh;
my %infa = ();
open my $ffh, '<', $fasta or die $!;
while(<$ffh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	$infa{$1} = 1;
    }
}
close $ffh;

my $ianf = 0;
print "IN ANNOTATION, NOT IN FASTA:\n";
foreach my $ks (sort keys %inann){
    unless(exists $infa{$ks}){
	print "\t".$ks."\n";
	$ianf += 1;
    }
}
my $ifna = 0;
print "IN FASTA, NOT IN ANNOTATION:\n";
foreach my $ks (sort keys %infa){
    unless(exists $inann{$ks}){
	print "\t".$ks."\n";
	$ifna += 1;
    }
}
print "In annotation, not in fasta:\t$ianf\nIn fasta, not in annotation:\t$ifna\n";
