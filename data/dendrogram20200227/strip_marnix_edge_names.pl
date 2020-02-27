#!/bin/env perl

use strict;
use warnings;

open my $cfh, '<', 'contig_border_table.txt' or die $!;
while(<$cfh>){
    chomp;
    my ($p, $left, $right) = split(/\t/, $_);
    next if($left eq '?' || $right eq '?');
    if($left==1 || $right==1){
	#my @paths = split(/\|/, $p);
	#foreach my $path (@paths){
	#    print "$path\n";
	#}
	print "$p\n";
    }
}
close $cfh;
