#!/bin/env perl

use strict;
use warnings;

foreach(@ARGV){
    my $afa = $_;
    $afa =~ s/\.phy$/\.afa/;
    open my $pfh, '<', $_ or die $!;
    open my $afh, '>', $afa or die $!;
    my $first = 1;
    while(<$pfh>){
	chomp;
	if($first){
	    $first = 0;
	    next;
	}else{
	    my ($header, $aln) = split(/\s+/, $_);
	    print $afh '>'.$header."\n$aln\n";
	}
    }
    close $pfh;
    close $afh;
}
