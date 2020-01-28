#!/bin/env perl

use strict;
use warnings;

my ($lo, $hi) = (undef, undef);
while(<>){
    chomp;
    my ($name, @score) = split(',', $_);
    foreach my $s (@score){
	$lo = $s unless(defined $lo);
	$hi = $s unless(defined $hi);
	$lo = $s if($s < $lo);
	$hi = $s if($s > $hi);
    }
}
print "LOW:\t$lo\nHIGH:\t$hi\n";
