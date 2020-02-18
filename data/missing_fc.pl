#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my %f = ();
open my $fc, '<', 'KS_fullclades.tsv' or die $!;
while(<$fc>){
    chomp;
    my ($ks, $clade) = split(/\t/, $_);
    $f{$ks} = $clade;
}
close $fc;

my $t = new Bio::SeqIO(-file=>'KS_rawseq_pred_training_transATPKS.txt', -format=>'fasta');
while(my $seq = $t->next_seq){
    print $seq->id."\n" unless(exists $f{$seq->id});
}
