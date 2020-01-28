#!/bin/env perl

use strict;
use warnings;

my @dorder = ();
foreach my $a (1..38){
    push @dorder, 'KS'.$a;
}
foreach my $b (1..4){
    push @dorder, 'KS'.$b.'b';
}

my %clade = ();
my $last = undef;
open my $dfh, '<', 'completeBGCs_mc.faa' or die $!;
while(<$dfh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	$last=$1;
    }else{
	@{$clade{$last}} = split(/,/, $_);
    }
}
close $dfh;

my %full = ();
open my $kfh, '<', 'completeBGCs_mc_index.faa' or die $!;
while(<$kfh>){
    chomp;
    if($_ =~ m/^>(.+)/){
	$last=$1;
    }else{
	my @k = split(/,/, $_);
	my $i = 0;
	foreach my $ks (@k){
	    $full{$last}{$ks} = ${$clade{$last}}[$i];
	    $i += 1;
	}
    }
}
close $kfh;

print "\t".join("\t", 'compd_KS_count', @dorder)."\n";
foreach my $pathway (sort keys %full){
    print join("\t", $pathway, scalar(keys %{$full{$pathway}}));
    foreach my $ks (@dorder){
	print "\t";
	if(exists $full{$pathway}{$ks}){
	    print $full{$pathway}{$ks};
	}else{
	    print 'NA';
	}
    }
    print "\n";
}
