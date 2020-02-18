#!/bin/env perl

use strict;
use warnings;

foreach my $faa (@ARGV){
    my %f = ();
    open my $ffh, '<', $faa or die $!;
    while(<$ffh>){
	chomp;
	if($_ =~ m/^>(.+)/){
	    $f{$1} += 1;
	}
    }
    close $ffh;
    foreach my $h (keys %f){
	if($f{$h} > 1){
	    print join("\t", $faa, $f{$h}, $h)."\n";
	}
    }
}
