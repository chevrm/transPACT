#!/bin/env perl

use strict;
use warnings;

my $i = 'annotateKS_per_pred_transATPKS_mc_manualderep.txt';
my $o = 'annotateKS_per_pred_transATPKS_mc_manualderep_20181029.txt';

my %blacklist = ();
open my $bfh, '<', 'annotation_dups_20181029.txt' or die $!;
while(<$bfh>){
    next if($_ =~ m/^\s/);
    chomp;
    $blacklist{$_} = 1;
}
close $bfh;

open my $ifh, '<', $i or die $!;
open my $ofh, '>', $o or die $!;
while(<$ifh>){
    chomp;
    my $ln = $_;
    $_ =~ s/^(\S+).+/$1/;
    unless(exists $blacklist{$_}){
	print $ofh "$ln\n";
    }
}
close $ifh;
close $ofh;

my %badlist = ();
open my $ofh2, '<', $o or die $!;
while(<$ofh2>){
    chomp;
    next if($_ =~ m/^\s/);
    $_ =~ s/^(\S+).+/$1/;
    
    my @q = split(/\|/, $_);
    unless($q[-1] =~ m/^c\d+$/){
	foreach my $acc (@q){
	    if($acc =~ m/^[A-Z].+\d$/){
		$badlist{$acc} = 1;
	    }
	}
    }    
}
close $ofh2;

my $j = 0;
my %whitelist = (
    'HM071004' => 1,
    'CP003468' => 1,
    'CP001614' => 1
    );
foreach my $bad (sort keys %badlist){
    next if(exists $whitelist{$bad});
    chomp(my $n = `grep -c $bad $o`);
    if($n > 1){
	chomp(my $t = `grep $bad $o`);
	my @ln = split(/\n/, $t);
	my ($q, $count, @rest) = split(/\s+/, $ln[0]);
	my $same = 1;
	foreach my $l (@ln){
	    my ($q2, $count2, @rest2) = split(/\s+/, $l);
	    if($count != $count2){
		$same = 0;
	    }
	}
	if($same == 1){
	    $j += 1;
	    print join("\t", $bad, "count=$n")."\n";
	    print join("\n", @ln)."\n";
	    print "\n\n";
	}
    }
}
print STDERR $j."\n";
