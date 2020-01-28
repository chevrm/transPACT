#!/bin/env perl

use strict;
use warnings;

open my $cfh, '<', 'KS_rawseq_pred_training_transATPKS.clades.tsv' or die $!;
my ($mult, $tot) = (0,0);
while(<$cfh>){
    chomp;
    my ($ks, $clade) = split("\t", $_);
    $mult += 1 if($clade =~ m/\|/);
    $tot += 1;
}
close $cfh;
print join("\t", $mult, $tot, $mult/$tot)."\n";
