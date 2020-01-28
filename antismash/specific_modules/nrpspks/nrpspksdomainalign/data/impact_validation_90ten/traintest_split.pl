#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use List::Util qw/shuffle/;

my $splits = 10;
my $trainsz = 583;
my $datadir = '/home/mchevrette/git/antismash-transat/antismash/specific_modules/nrpspks/nrpspksdomainalign/data/';
my $annotation = $datadir.'annotation_mc_20180104.tsv';
my $fasta = $datadir.'647KS_mcformat.afa';

## Load seqs
my %fa = ();
my $f = new Bio::SeqIO(-format=>'fasta', -file=>$fasta);
while(my $seq=$f->next_seq){
    my $nogap = $seq->seq;
    $nogap =~ s/-//g;
    $fa{$seq->id} = $nogap;
}

## Name check
my %a = ();
open my $afh, '<', $annotation or die $!;
while(<$afh>){
    chomp;
    next if($_ =~ m/^#/);
    my($leaf, $clade, $desc) = split(/\t/, $_);
    $a{$leaf} = 1;
}
close $afh;
my $good = 0;
foreach my $sid (keys %fa){
    if(exists $a{$sid}){
	$good+=1;
    }else{
	die STDERR $sid." not found.\n";
    }
}
#print STDERR $good." sequences correct.\n";

for(my $s = 1;$s<=$splits;$s+=1){
    my $splitname = 'traintest_split_';
    if($s < 10){
	$splitname .= '0'.$s;
    }else{
	$splitname .= $s;
    }
    system("rm -r $splitname") if(-d $splitname);
    mkdir $splitname;
    my $trainfile = $splitname.'/'.$splitname.'_train.faa';
    my $testfile = $splitname.'/'.$splitname.'_test.faa';
    my @names = shuffle(keys %fa);
    my @trainarr = @names[0..$trainsz-1];
    my @testarr = @names[$trainsz..$#names];
    open my $train, '>', $trainfile or die $!;
    foreach my $seq (@trainarr){
	print $train '>'.$seq."\n".$fa{$seq}."\n";
    }
    close $train;
    open my $test, '>', $testfile or die $!;
    foreach my $seq (@testarr){
	print $test '>'.$seq."\n".$fa{$seq}."\n";
    }
    close $test;
}
