#!/bin/env perl

use strict;
use warnings;

my %c = ();
while(<>){
    chomp;
    if($_ =~ m/^>(.+)/){
	my @h = split(/\|/, $1);
	my $name = join('|', @h[0..1]);
	my $ks = $h[-1];
	$c{$name}{$ks} += 1;
    }
}

foreach my $name (keys %c){
    foreach my $ks (keys %{$c{$name}}){
	print join("\t", $c{$name}{$ks}, $name, $ks)."\n";
    }
}
