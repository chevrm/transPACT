#!/bin/env perl

use strict;
use warnings;

my %f = ();
while(<>){
    chomp;
    my ($ks, $clade) = split(/\t/, $_);
    if(exists $f{$ks} && $clade ne $f{$ks}){
	print join("\t", $ks, $f{$ks}, $clade)."\n";
    }
    $f{$ks} = $clade;
}
foreach my $ks (sort keys %f){
    #print join("\t", $ks, $f{$ks})."\n" if($f{$ks} > 1);
}
