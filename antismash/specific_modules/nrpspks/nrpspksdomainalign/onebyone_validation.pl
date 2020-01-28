#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;


foreach my $v (1..10){
    $v = '0'.$v if($v < 10);
    my $d = './data/impact_validation/traintest_split_'.$v.'/';
    my $testfaa = $d.'traintest_split_'.$v.'_test.faa';
    my $traintree = $d.'RAxML_bestTree.traintest_split_'.$v.'_train_hmmalign_raxml.tre';
    my $trainafa = $d.'traintest_split_'.$v.'_train_hmmalign.afa';
    my $out = $d.'KS_cladevalidation_split_'.$v.'.tsv';
    my $faa = new Bio::SeqIO(-file=>$testfaa, -format=>'fasta');
    chomp(my $tot = `grep -c '>' $testfaa`);
    my $n = 1;
    system("rm $out") if(-e $out);
    while(my $seq = $faa->next_seq){
	open my $tfh, '>', 'tmp.faa' or die $!;
	print $tfh '>'.$seq->id."\n".$seq->seq."\n";
	close $tfh;
	system("python2 substrate_from_faa_validation.py tmp.faa $traintree $trainafa >> $out");
	print join("\t", 'split'.$v, $n, $tot, sprintf("%.3f", $n/$tot), $seq->id)."\n";
	$n += 1;
    }
}


