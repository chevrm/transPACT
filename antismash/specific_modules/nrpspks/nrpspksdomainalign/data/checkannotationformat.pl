#!/bin/env perl

use strict;
use warnings;

my %seen = ();
while(<>){
    chomp;
    next if($_ =~ m/^#/);
    my ($name, @rest) = split(/\t/, $_);
    $seen{$name} += 1;
}
foreach my $n (sort keys %seen){
    print $seen{$n}."\t$n\n" unless($seen{$n} == 2);
}
