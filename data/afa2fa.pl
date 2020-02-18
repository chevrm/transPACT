#!/bin/env perl

use strict;
use warnings;

while(<>){
    chomp;
    unless($_ =~ m/^>/){
	$_ =~ s/-//g;
    }
    print $_."\n";
}
