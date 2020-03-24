#!/bin/env perl

use strict;
use warnings;


my %kscount = ();
my %clustks = ();
open my $cfh, '<', 'KS_precomputed_1405_hmmalign_trimmed_renamed.clades.tsv' or die $!;
my $n = 0;
while(<$cfh>){
    chomp;
    my ($query, $clade, $ann) = split(/\t/, $_);
    $query =~ s/^q_//;
    my @q = split(/\|/, $query);
    my $ks = $q[-1];
    $ks =~ s/-$//;
    $ks = 'KS5' if($ks eq 'KS55'); # 9-methylstreptimidone_KS55
    my $clust = join('|', @q[0..$#q-1]);
    $ks = 'KS17' if($ks eq 'KS' && $clust eq 'phormidolide');
    if($query =~ m/OUTGROUP/){
	$ks = 'KS1';
	if($query =~ m/2qo3/){
	    $clust = 'OUTGROUP_2qo3|chainA_EryKS3';
	}elsif($query =~ m/2hg4/){
	    $clust = 'OUTGROUP_2hg4|chainA_EryKS5';
	}
    }
    if($clade =~ m/^Clade/){
	$kscount{$clust} += 1;
    }
    $clustks{$clust}{$ks} = $clade;
    $n += 1;
}
close $cfh;
## Reg = 1..38
## b = 1..4
my @order = ();
foreach my $n (1..38){
    push @order, 'KS'.$n;
}
foreach my $n (1..4){
    push @order, 'KS'.$n.'b';
}
print "\t".join("\t", @order)."\n";
foreach my $clust (sort keys %kscount){
    next unless($kscount{$clust} > 1 || $clust =~ m/^OUTGROUP/);
    my $ksstring = '';
    my $tot = scalar(keys %{$clustks{$clust}});
    foreach my $ks (@order){
	if(exists $clustks{$clust}{$ks}){
	    if($clustks{$clust}{$ks} =~ m/^Clade/ || $clustks{$clust}{$ks} eq 'clade_not_conserved'){
		if($ksstring eq ''){
		    $ksstring = $clustks{$clust}{$ks};
		}else{
		    $ksstring .= "\t".$clustks{$clust}{$ks};
		}
	    }#else{
		#last;
	    #}	
	}
    }
    foreach my $n ($tot..scalar(@order)-1){
	$ksstring .= "\t".'NA';
    }
    my ($o, $g) = (scalar(@order), scalar(split(/\t/, $ksstring)));
    die "You dont know how to count\n$clust\n$ksstring\n$o\t$g\n" unless($o == $g);
    print "$clust\t$ksstring\n";
}
