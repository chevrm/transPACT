#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my %d = ();
my %col = ();
my %seen = ();
my %c = ('Other' => '#d9d9d9',
	 'bOH' => '#a6cee3',
	 'eDB' => '#1f78b4',
	 'double_bonds' => '#b2df8a',
	 'Starter' => '#33a02c',
	 'non_elongating' => '#fb9a99',
	 'amino_acids' => '#e31a1c',
	 'furan|pyran' => '#fdbf6f',
	 'exometh|red_nMe' => '#ff7f00',
	 'zDB' => '#cab2d6'
);
open my $afh, '<', '../dendrogram20190514/Annotation_MC_final.txt' or die $!;
while(<$afh>){
    chomp;
    next if($_ =~ m/^#/);
    my ($leaf, $clade, $desc) = split(/\t/, $_);
    $clade =~ s/Clade_//;
    $clade = '0'.$clade if($clade < 10);
    $desc =~ s/^\s*//;
    $desc =~ s/\s*$//;
    $desc =~ s/furan\s/furan/;
    $desc =~ s/amino\s/amino_/;
    if($desc =~ m/double.*bonds/i){
	$desc = 'double_bonds';
    }elsif($desc =~ m/b.*OH/i){
	$desc = 'bOH';
    }elsif($desc =~ m/e.*DB/i){
	$desc = 'eDB';
    }elsif($desc =~ m/Starter/i){
	$desc = 'Starter';
    }elsif($desc =~ m/non_elongating/){
	$desc = 'non_elongating';
    }
    $desc = 'Other' unless(exists $c{$desc});
    $d{$clade} = $c{$desc};
}
close $afh;
$d{'NC'} = $c{'Other'};
$d{'NA'} = $c{'Other'};

my %p = ();
open my $bfh, '<', '../dendrogram20190514/KS_precomputed_1405_hmmalign_trimmed_renamed.annotationtable.tsv';
while(<$bfh>){
    chomp;
    next if($_ =~ m/^\s/);
    my ($pathway, @clades) = split(/\t/, $_);
    my @a = ();
    foreach my $c (@clades){
	next if($c eq 'NA');
	if($c eq 'clade_not_conserved'){
	    $c = 'NC';
	}else{
	    $c =~ s/Clade_//;
	    $c = '0'.$c if($c < 10);
	}
	push @a, $c;
    }
    $p{$pathway} = join('-', @a);
}
close $bfh;

my %isce = ();
open my $cfh, '<', 'leaftbl.tsv' or die $!;
while(<$cfh>){
    chomp;
    my ($p, $ce) = split(/\t/, $_);
    if($ce eq 'Yes'){
	$isce{$p} = 1;
    }
}
close $cfh;


print join("\n", 'DATASET_DOMAINS', 'SEPARATOR COMMA', 'DATASET_LABEL,ks', "COLOR,#ff0000\n", 'DATA')."\n";
print STDERR join("\n", 'DATASET_BINARY', 'SEPARATOR COMMA', 'DATASET_LABEL,Contig edge', 'COLOR,#ff0000', 'FIELD_SHAPES,1', 'FIELD_LABELS,f1', 'DATA')."\n";
my $sz = 100;
my $t = new Bio::TreeIO(-file=>'./upgma_trim_species.nwk', -format=>'newick');
while(my $tree = $t->next_tree){
    foreach my $leaf ($tree->get_leaf_nodes){
	next if($leaf->id =~ m/^Inner/);
	my @l = split(/\|/, $leaf->id);
	my $last2 = $l[-1];
	if(defined $l[-2]){
	    $last2 = join('|', $l[-2], $l[-1]);
	}
	die join("\n", $leaf->id, $last2)."\n" unless(exists $p{$last2});
	my @ks = split('-', $p{$last2});
	#die $leaf->id."\n" if (scalar(@ks)==0);
	my $ln = $leaf->id.','.scalar(@ks)*$sz.',';
	my $start = 0;
	foreach my $k (@ks){
	    die "DIED: $k\n" unless(exists $d{$k});
	    $ln .= join('|', 'EL', $start, $start+$sz, $d{$k}, $k).',';
	    $start += $sz
	}
	$ln =~ s/,$//;
	print "$ln\n";
	if(exists $isce{$last2}){
	    print STDERR join(',', $leaf->id,1)."\n";
	}else{
	    print STDERR join(',', $leaf->id,-1)."\n";
	}
    }
}
