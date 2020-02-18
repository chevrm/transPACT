#!/bin/env perl

use strict;
use warnings;

my %seen = ();
while(<>){
    chomp;
    my @l = split(/,/, $_);
    $seen{$l[2]} += 1;
}

foreach my $l (sort {$seen{$b} <=> $seen{$a}} keys %seen){
    print join("\t", $seen{$l}, $l)."\n";
}
