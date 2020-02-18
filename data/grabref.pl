#!/bin/env perl

use strict;
use warnings;

my %r = ();
my @lines = ();
while(<>){
    chomp;
    push @lines, $_;
    my @l = split(/\t/, $_);
    if($l[0] =~ m/^.+\|referencetree$/){
	my ($f, @rest) = split(/\|/, $l[0]);
	$r{$f} = 1;
    }
}

foreach my $l (@lines){
    my ($h, @rest) = split(/\t/, $l);
    my @f = split(/\|/, $h);
    if(exists $r{$f[0]}){
	print $l."\n";
    }
}
