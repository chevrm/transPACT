#!/bin/env perl

use strict;
use warnings;

my @j = (0.22, 0.00);
my @d = (0.76, 0.32);
my @g = (0.00, 0.00);
my @a = (0.02, 0.68);

for(my $i=0;$i<scalar(@j);$i += 1){
    my $dir = 'j'.$j[$i].'g'.$g[$i].'d'.$d[$i].'a'.$a[$i];
    $dir =~ s/0\.//g;
    print $dir."\n";
    if(-d $dir){
	next;
    }else{
	mkdir $dir;
	system("python generate_dendrogram_userweights.py $j[$i] $g[$i] $d[$i] $a[$i]");
	system("perl tree_rn_species.pl");
	system("perl generate_itol.pl > itol.txt");
	system("mv *txt $dir");
	system("mv *nwk $dir");
    }
}
