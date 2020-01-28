#!/bin/env perl

use strict;
use warnings;

my %h = ();
while(<>){
    chomp;
    if($_ =~ m/^>(.+)\s*/){
	$h{$1} += 1;
    }
}
foreach my $head (sort keys %h){
    if($h{$head} > 1){
	print join("\t", $h{$head}, $head)."\n";
    }
}
