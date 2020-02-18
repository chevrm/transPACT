#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $out = 'KS_precomputed_hmmalign.afa';
system("rm $out") if(-e $out);
my $hmm = 'KS_rawseq_pred_training_transATPKS.hmm';
my $faa = new Bio::SeqIO(-format=>'fasta', -file=>'KS_rawseq_pred_training_transATPKS.txt');
while(my $seq = $faa->next_seq){
    open my $tfh, '>', 'tmp.faa' or die $!;
    print $tfh '>'.$seq->id."\n".$seq->seq."\n";
    close $tfh;
    system("hmmalign $hmm tmp.faa > tmp.stk");
    system("python stk2aln.py >> $out");
}
system("rm tmp.faa tmp.stk");
