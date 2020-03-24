#!/bin/env perl

use strict;
use warnings;

my @j = (0.00);
my @d = (0.32);
my @a = (0.68);

for(my $i=0;$i<scalar(@j);$i += 1){
    my $dir = 'j'.$j[$i].'d'.$d[$i].'a'.$a[$i];
    $dir =~ s/0\.//g;
    print $dir."\n";
    mkdir $dir unless(-d $dir);
    system("python generate_dendrogram_userweights.py $j[$i] $d[$i] $a[$i]") unless(-e 'upgma.nwk');
    system("perl tree_rn_species.pl");
    system("perl generate_itol.pl > itol.txt");
    system("mv *txt $dir");
    system("mv *nwk $dir");
}
