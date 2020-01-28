#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $fa = new Bio::SeqIO(-file=>shift, -format=>'fasta');

my $s1 = $fa->next_seq;
my $s2 = $fa->next_seq;

my @seq1 = split(//, $s1->seq);
my @seq2 = split(//, $s2->seq);

my ($gap, $match) = (0,0);
for(my $i=0;$i<scalar(@seq1);$i++){
    if($seq1[$i] eq '-'){
	$gap++;
    }elsif($seq2[$i] eq '-'){
	$gap++;
    }elsif($seq1[$i] eq $seq2[$i]){
	$match++;
    }
}
my $tot = scalar(@seq1) - $gap;
print join("\t", $match, $tot, $match/$tot)."\n";
