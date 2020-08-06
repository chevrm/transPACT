#!/bin/env perl

use strict;
use warnings;

my %s = ();
my %d = ();
open my $sfh, '<', 'strands.tsv' or die $!;
while(<$sfh>){
    chomp;
    my ($bgc, $ks, $start, $end, $strand) = split(/\t/, $_);
    push @{$s{$bgc}}, $ks;
    $d{$bgc} += $strand;
}
close $sfh;

my %cladeof = ();
open my $pfh, '<', 'pks_pkslike_ks.trim.tsv' or die $!;
while(<$pfh>){
    chomp;
    my ($longks, $clade, $ann) = split(/\t/, $_);
    my ($bgc, $ks) = ($longks, $longks);
    $bgc =~ s/(.+k\d{3}).+/$1/;
    $ks =~ s/.+k\d{3}_(.+)/$1/;
    $cladeof{$bgc}{$ks} = $clade;
}
close $pfh;

## Add precomputed outgroups
$cladeof{'OUTGROUP_2hg4|chainA_EryKS5'}{'KS1'} = 'Clade_40';
$d{'OUTGROUP_2hg4|chainA_EryKS5'} = 1;
@{$s{'OUTGROUP_2hg4|chainA_EryKS5'}} = ('KS1');
$cladeof{'OUTGROUP_2qo3|chainA_EryKS3'}{'KS1'} = 'Clade_40';
$d{'OUTGROUP_2qo3|chainA_EryKS3'} = 1;
@{$s{'OUTGROUP_2qo3|chainA_EryKS3'}} = ('KS1');

## Reg = 1..38
## b = 1..4
my @order = ();
foreach my $n (1..38){
    push @order, 'KS'.$n;
}
print "\t".join("\t", @order)."\n";
foreach my $bgc (sort keys %cladeof){
    my $n = scalar(@{$s{$bgc}});
    next unless($n > 1 || $bgc =~ m/OUTGROUP/);
    print $bgc;
    my $l = 38;
    if($d{$bgc} >= 0){
	foreach my $ks (@{$s{$bgc}}){
	    print "\t".$cladeof{$bgc}{$ks};
	    $l = $l-1;
	}
    }else{
	foreach my $ks (reverse @{$s{$bgc}}){
	    print "\t".$cladeof{$bgc}{$ks};
	    $l = $l-1;
	}
    }
    while($l > 1){
	print "\tNA";
	$l = $l-1;
    }
    print "\n";
}
