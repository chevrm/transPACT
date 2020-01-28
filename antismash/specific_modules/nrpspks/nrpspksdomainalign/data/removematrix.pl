#!/bin/env perl

use strict;
use warnings;

my $from = 'annotateKS_per_pred_transATPKS_mc.txt';
my $to = 'annotateKS_per_pred_transATPKS_mc_manualderep.txt';
my $rm = 'toremovefrommatrix.txt';

my %blacklist = ();
open my $rfh, '<', $rm or die $!;
while(<$rfh>){
    chomp;
    $blacklist{$_} = 1;
}
close $rfh;

open my $ifh, '<', $from or die $!;
open my $ofh, '>', $to or die $!;
while(<$ifh>){
    chomp;
    my ($h, @rest) = split(/\t/, $_);
    print $ofh $_."\n" unless(exists $blacklist{$h});
}
close $ifh;
close $ofh;
