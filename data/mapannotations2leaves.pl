#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $ann_file = 'annotation_eh_20171106.txt';
my $tree_file = 'RAxML_bestTree.647KS_RAxML.tre';

my %leaf = ();
my $t = new Bio::TreeIO(-file=>$tree_file, -format=>'newick');
while(my $tree = $t->next_tree){
    foreach my $l ($tree->get_leaf_nodes){
	$leaf{$l->id} = 1;
    }
}

#print "LEAVES FOUND:\t".scalar(keys %leaf)."\n";

my %a = ();
open my $afh, '<', $ann_file or die $!;
my @noclade = ();
while(<$afh>){
    next if($_ =~ m/^\s+$/);
    chomp;
    my ($name, $clade, $ann) = ($_, $_, $_);
    $name =~ s/^\s*(.+)\s+Clade_.+$/$1/;
    $name =~ s/\s+$//;
    $name =~ s/\s+/_/g;
    if($clade =~ m/^.+(Clade_\d+)\s+(.+)$/){
	$clade = $1;
	$ann = $2;
	$ann =~ s/\s+$//;
	$ann =~ s/\//|/g;
	$ann =~ s/\s+/_/g;
    }else{
	$clade = 'no_clade';
	$ann = 'no_clade';
	push @noclade, $name;
    }
    $a{$name}{'c'} = $clade;
    $a{$name}{'a'} = $ann;
}
close $afh;

print join("\t", '#Leaf_name', 'Clade_name', 'Clade_desc')."\n";
my (@notintree,@notinann) = ((),());
foreach my $n (sort keys %a){
    print join("\t", $n, $a{$n}{'c'}, $a{$n}{'a'})."\n";
    unless(exists $leaf{$n}){
	push @notintree, $n;
    }
}
foreach my $l (sort keys %leaf){
    unless(exists $a{$l}){
	push @notinann, $l;
    }
}

#print STDERR "NOT IN TREE:\t".scalar(@notintree)."\nNOT IN ANN:\t".scalar(@notinann)."\n";
print join("\n", sort @notinann)."\n";
#print join("\n", sort @noclade)."\n";
