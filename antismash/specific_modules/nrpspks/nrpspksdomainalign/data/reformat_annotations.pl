#!/bin/env perl

use strict;
use warnings;

print join("\t", '#Leaf_name', 'Clade_name', 'Clade_desc')."\n";
my $sc = 1;
while(<>){
    my $t = $_;
    chomp;
    $_ =~ s/\r//g;
    next if($_ =~ m/^\s*$/);
    my ($leaf_name, $clade_name, $clade_desc) = split(/\t/, $_);
    $leaf_name =~ s/\s+/_/g;
    if(defined $clade_desc && $clade_name =~ m/^Clade/){
	my @d = split(/\//, $clade_desc);
	my $new_desc = join('|', sort {$a cmp $b} @d);
	print join("\t", $leaf_name, $clade_name, $new_desc)."\n";
    }else{
	my $new_desc = $leaf_name;
	$new_desc =~ s/^.+KS\d+_(.+)$/$1/;
	print join("\t", $leaf_name, 'AutoClade_'.$sc, $new_desc)."\n";
	$sc += 1;
    }

}
