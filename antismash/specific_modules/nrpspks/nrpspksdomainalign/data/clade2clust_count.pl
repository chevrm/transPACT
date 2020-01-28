#!/bin/env perl

use strict;
use warnings;

my $clade = './KS_fullclades.tsv';

my %c2k = ();

open my $ifh, '<', $clade or die $!;
while(<$ifh>){
    chomp;
    my ($name, $clade) = split(/\t/, $_);
    $clade = 'clade_not_conserved' if($clade =~ m/\|/);
    my @n = split(/\|/, $name);
    my $k = $n[0];
    if(scalar(@n) > 1 && $n[1] ne ''){
	$k .= '|'.$n[1];
    }
    my $ksnum = $n[-1];
    $ksnum = $n[-2] if($ksnum eq 'KS');
    $c2k{$k}{$ksnum} = $clade;
}
close $ifh;

my ($u, $p) = (0, 0);
foreach my $k (sort keys %c2k){
    next if($k =~ m/OUTGROUP/);
    my $unordered = 0;
    my $lastks = 0;
    my (@ksa, @c) = ((), ());;
    foreach my $ks (sort kssort keys %{$c2k{$k}}){
	my $ksnum = $ks;
	$ksnum =~ s/^\D+(\d+)\D*$/$1/;
	unless($ksnum == $lastks || $ksnum == $lastks+1){
	    $unordered = 1;
	    $u += 1;
	}
	$lastks = $ksnum;
	push @ksa, $ks;
	push @c, $c2k{$k}{$ks};
    }
    if($unordered || scalar(@ksa) <= 1){
	$p += 1;
    }
    #unless($unordered || scalar(@ksa) <= 1){ ## UNCOMMENT IF YOU WANT ORDERED ONLY!
	#print '>'.$k."\n".join(',', @c)."\n";
	#print STDERR '>'.$k."\n".join(',', @ksa)."\n";
    #}
}

print "$u unordered, $p pathways.\n";

sub kssort{
    my ($ks1, $ks2) = ($a, $b);
    if($ks1 =~ m/^KS/ && $ks2 =~ m/^KS/){
	$ks1 =~ s/^KS//;
	$ks2 =~ s/^KS//;
	my $ks1n = $ks1;
	$ks1n =~ s/^(\d+)\D+/$1/;
	my $ks2n = $ks2;
	$ks2n =~ s/^(\d+)\D+/$1/;
	if($ks1n != $ks2n){
	    return $ks1n <=> $ks2n;
	}else{
	    my ($ks1l, $ks2l) = ('a', 'a');
	    if($ks1 =~ m/^\d+(\S)$/){
		$ks1l = $1;
	    }
	    if($ks2 =~ m/^\d+(\S)$/){
		$ks2l = $1;
	    }
	    return $ks1l cmp $ks2l;
	}
    }else{
	return $ks1 cmp $ks2;
    }
}
