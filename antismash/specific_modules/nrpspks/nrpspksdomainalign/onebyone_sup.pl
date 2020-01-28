#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $faf = './data/647KS_mhm_2.fasta';
my $faa = new Bio::SeqIO(-file=>$faf, -format=>'fasta');
chomp(my $tot = `grep -c '>' $faf`);
my $n = 1;
my $out = './data/KS_fullclades_sup.tsv';
system("rm $out") if(-e $out);
while(my $seq = $faa->next_seq){
    print $seq->id."\n";
    open my $tfh, '>', 'tmp.faa' or die $!;
    print $tfh '>'.$seq->id."\n".$seq->seq."\n";
    close $tfh;
    system("python substrate_from_faa.py tmp.faa >> $out");
    print "$n of $tot complete\n";# if($n % 500 == 0);
    $n += 1;
}
