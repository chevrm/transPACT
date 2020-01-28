#!/bin/env perl

use strict;
use warnings;

my @header = ('Leafname', 'GenbkID', 'KS_index', 'ClusterID', 'Clade');
print "\t".join("\t", @header)."\n";
my $n = 0;
my %cladeof = ();
open my $cfh, '<', 'KS_fullclades.tsv' or die $!;
while(<$cfh>){
    chomp;
    my ($ks, $clade) = split(/\t/, $_);
    $cladeof{$ks} = $clade;
}
close $cfh;

foreach my $ks (sort keys %cladeof){
    my @h = split(/\|/, $ks);
    my $gid = $h[0];
    my $ksi = $h[-1];
    my $ci = $h[1];
    print STDERR join("\t", $ks, $ci)."\n";
    print join("\t", $n, $ks, $gid, $ksi, $ci, $cladeof{$ks})."\n";
    $n+=1;
}
