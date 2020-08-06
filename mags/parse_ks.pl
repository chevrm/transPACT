#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

foreach my $gbf (@ARGV){
    my @p = split(/\//, $gbf);
    my $clust = $p[-1];
    $clust =~ s/region(\d+).+/k$1/;
    my $gb = new Bio::SeqIO(-file=>$gbf, -format=>'genbank');
    while(my $seq = $gb->next_seq){
	#print $seq->desc."\n";
	for my $feat ($seq->get_SeqFeatures){
	    if ($feat->primary_tag eq 'aSDomain'){
		my ($id, $pseq, $type) = (undef, undef, undef);
		foreach my $tag ($feat->get_tag_values('aSDomain') ){
		    $type = $tag;
		}
		next unless($type =~ m/^PKS_KS/ || $type =~ m/^mod_KS/ || $type =~ m/^hyb_KS/ || $type =~ m/^tra_KS/);
		foreach my $tag ($feat->get_tag_values('translation') ){
		    $pseq = $tag;
		}
		foreach my $tag ($feat->get_tag_values('label') ){
		    $id = $tag;
		}
		print '>'.join('_', $clust, $id)."\n$pseq\n";
	    }
	}
    }
}
