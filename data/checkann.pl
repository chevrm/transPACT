#!/bin/env perl

use strict;
use warnings;

my %l = ();
while(<>){
    chomp;
    next if($_ =~ m/^#/);
    my ($leaf, @rest) = split(/\t/, $_);
    $l{$leaf} = 1;
}
open my $a, '<', '647KS_mcformat.afa' or die $!;
while(<$a>){
    chomp;
    if($_ =~ m/^>(.+)$/){
	print $1."\n" unless(exists $l{$1});
    }
}
close $a;
