#!/bin/env perl

use strict;
use warnings;

my %seen = ();
while(<>){
    chomp;
    my @l = split(/,/, $_);
    if($l[-1] eq '?'){
	$seen{$l[2]} += 1;
    }
}

foreach my $l (sort {$seen{$b} <=> $seen{$a}} keys %seen){
    print join("\t", $l)."\n";
}
