#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my %l = ();
my $faa = new Bio::SeqIO(-file=>shift, -format=>'fasta');
while(my $seq = $faa->next_seq){
    $l{$seq->length} += 1;
}

print join("\t", 'Length', 'Count')."\n";
foreach my $len(sort keys %l){
    print join("\t", $len, $l{$len})."\n";
}
