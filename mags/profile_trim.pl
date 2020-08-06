#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $ref = '../data/647KS_mcformat.afa';
chomp(my $refn = `grep -c '>' $ref`);
my $qf = shift;
my $nf = $qf;
$nf =~ s/\.faa$/.trim.faa/;

open my $nfh, '>', $nf or die $!;
my $fa = new Bio::SeqIO(-file=>$qf, -format=>'fasta');
while(my $seq = $fa->next_seq){
    open my $tfh, '>', 'tmp.faa' or die $!;
    print $tfh '>'.$seq->id."\n".$seq->seq."\n";
    close $tfh;
    system("muscle -profile -in1 $ref -in2 tmp.faa -out tmp.afa");
    my %gh = ();
    my %gt = ();
    my $afa = new Bio::SeqIO(-file=>'tmp.afa', -format=>'fasta');
    while(my $aseq = $afa->next_seq){
	next if($aseq->id eq $seq->id);
	foreach my $n (1..10){
	    $gh{$n} += 1 if($aseq->seq =~ m/^-{$n}/);
	    $gt{$n} += 1 if($aseq->seq =~ m/-{$n}$/);
	}
    }
    my ($f, $t) = (0,0);
    foreach my $n (1..10){
	$f = $n if(exists $gh{$n} && $gh{$n} == $refn);
	$t = $n if(exists $gt{$n} && $gt{$n} == $refn);
    }
    my $trimseq = $seq->seq;
    $trimseq =~ s/^\w{$f}(\w+)\w{$t}$/$1/;
    print $nfh '>'.$seq->id."\n".$trimseq."\n";
}
close $nf;
