#!/bin/env perl

use strict;
use warnings;

my %f = ();
open my $ffh, '<', shift or die $!;
while(<$ffh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	$f{$1} = 1;
    }
}
close $ffh;

open my $sfh, '<', shift or die $!;
while(<$sfh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	unless(exists $f{$1}){
	    print $1."\n";
	}
    }
}
close $sfh;
