#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $reftree = 'RAxML_bestTree.649KS_sequences_hmmalign_raxml.tre';
my $ann = 'Annotation_Marnix_final.txt';

my %leaf = ();
my $in = new Bio::TreeIO(-file=>$reftree, -format=>'newick');
while(my $tree = $in->next_tree){
    foreach my $l ($tree->get_leaf_nodes){
        $leaf{$l->id} = 1;
    }
}

print join("\t", '#Leaf_name','Clade_name','Clade_desc')."\n";
open my $afh, '<', $ann or die $!;
while(<$afh>){
    chomp;
    my @ln = split(/\s+/, $_);
    next unless(scalar(@ln) > 1);
    my $n = 0;
    my ($k, $c) = (0,0);
    foreach my $e (@ln){
	if($e =~ m/^KS\d/){
	    $k = $n;
	}elsif($e =~ m/^Clade_\d/){
	    $c = $n;
	}
	$n += 1;
    }
    my $name1 = '';
    my $name2 = join('', @ln[$k+1..$c-1]);
    if($ln[$k-1] =~ m/^\d+$/){
	$name1 = join('_', @ln[0..$k-2]).$ln[$k-1].'_'.$ln[$k];
    }else{
	$name1 = join('_', @ln[0..$k]); 
    }
    if($ln[0] eq 'phormidolide'){
	$name1 = join('_', @ln[0..$k]); 
    }

    $name1 =~ s/^9_methylstreptimidone/9-methylstreptimidone/;
    $name1 =~ s/^Calyculin/calyculin/;
    $name1 =~ s/CalC_(\d+)/CalC$1/;
    $name1 =~ s/^methoxycarbonyl_myxopyronin/myxopyronin/;
    $name1 =~ s/^vinylogous_//;

    $name2 = '0bimod_aMeeDB' if($name2 eq '0bimodaMeeDB');
    $name2 = '0biomod_aMeeDB' if($name2 eq '0biomodaMeeDB');
    $name2 = 'eDB_zDB' if($name2 eq 'eDBzDB');
    $name2 = 'eDB_DA' if($name2 eq 'eDBDA');
    $name2 = 'DB_DA' if($name2 eq 'DBDA');
    $name2 = 'zDB_eDB' if($name2 eq 'zDBeDB');
    $name2 = 'zD_zDB' if($name2 eq 'zDzDB');
    $name2 = 'zDB_zDB' if($name2 eq 'zDBzDB');
    $name2 = 'bimod_eDB' if($name2 eq 'bimodeDB');    
    $name2 = '0Hbimod_bOH' if($name2 eq '0HbimodbOH');
    $name2 = '0bimod_bDOH' if($name2 eq '0bimodbDOH');
    $name2 = '0bimod_bOH' if($name2 eq '0bimodbOH');
    $name2 = 'Obimod_bOH' if($name2 eq 'ObimodbOH');
    $name2 = 'OHbimod_OH' if($name2 eq 'OHbimodOH');
    $name2 = 'obimod_bOH' if($name2 eq 'obimodbOH');
    $name2 = 'bimod_bOH' if($name2 eq 'bimodbOH');
    $name2 = 'bimod_bDOH' if($name2 eq 'bimodbDOH');
    $name2 = '0bimod_eDB' if($name2 eq '0bimodeDB');
    $name2 = 'biMod_eDB' if($name2 eq 'biModeDB');
    $name2 = 'bimod_0bLOH' if($name2 eq 'bimod0bLOH');
    $name2 = '0_bLOH_aLMebLOH' if($name2 eq '0bLOHaLMebLOH');
    $name2 = '0bimod_bDOMe' if($name2 eq '0bimodbDOMe');
    $name2 = '0bimod_aDMebDOH' if($name2 eq '0bimodaDMebDOH');
    $name2 = 'obimod_bLOH' if($name2 eq 'obimodbLOH');    

    
    $name2 = 'unusualStarters_Succ' if($name2 eq 'unusualStartersSucc');
    $name2 = 'unusualStarter_AMT' if($name2 eq 'unusualStarterAMT');
    $name2 = 'unusualstarter_AMT' if($name2 eq 'unusualstarterAMT');
    $name2 = 'aromaticStarter_AL' if($name2 eq 'aromaticStarterAL');
    $name2 = 'unusualStarters_acryloyl' if($name2 eq 'unusualStartersacryloyl');
    $name2 = 'unusualStarter_lactate' if($name2 eq 'unusualStarterlactate');
    $name2 = 'unusualstarter_lactate' if($name2 eq 'unusualstarterlactate');
    $name2 = '2biMod_OH' if($name2 eq '2biModOH');
    $name2 = 'bDOH_eDB' if($name2 eq 'bDOHeDB');
    
    my $name = join('_', $name1, $name2);
    $name =~ s/_$//;
    $name = 'malleilactone_burkholderic_acid_BurA1_KS1_unusualStarter' if($name eq 'malleilactone_burkholderic_acid_BurA1_KS2_unusualStarter');
    $name = 'tartrolon_TrtD_KS33_bimod_bOH' if($name eq 'tartrolon_TrtD_KS3_3bimodbOH');
    $name = 'basiliskamides_P615_14890_KS5_eDB' if($name eq 'basiliskamides_P61514890_KS5_eDB');
    $name = 'phormidolide_phmE_3_KS3_bOH' if($name eq 'phormidolide_phmE_3_KS3_bDOH');
    $name = 'corallopyronin_CorI1_KS1_unusualStarter_methoxycarbonyl' if($name eq 'corallopyronin_CorI1_KS1_unusualStarter');
    $name = 'calyculin_CalB1_KS3_KS0_bLOH' if($name eq 'calyculin_CalB_1_KS3_KS0_bLOH');
    $name = 'phormidolide_phmH_1_KS9_0aMebOH' if($name eq 'phormidolide_phmH_1_KS9_0aMeDOH');
    $name = 'dorrigocin_migrastatin_MgsE3_KS3_vinylogous' if($name eq 'dorrigocin_migrastatin_MgsE3_KS3');
    $name = 'bongkrekic_acid_BonA_KS11_GNATstarter' if($name eq 'bongkrekic_acid_BonA_KS1_1GNATstarter');
    $name = '2qo3_chainA_EryKS3_OUTGROUP' if($name eq '2qo3_chainAEryKS3OUTGROUP');
    $name = '2hg4_chainA_EryKS5_OUTGROUP' if($name eq '2hg4_chainAEryKS5OUTGROUP');
    $name = 'rhizopodin_RizB7_KS7_0bimodbDOH' if($name eq 'rhizopodin_RizB7_KS7_0bimod_bDOH');
    $name = 'bacillaene_Bamy_BaeL2_KS5_eD' if($name eq 'bacillaene_Bamy_BaeL2_KS5');
    $name = 'macrolactin_MlnG1_KS10_dH' if($name eq 'macrolactin_MlnG');
    $name = 'griseoviridin_SgvE22_KS4__eDB' if($name eq 'griseoviridin_SgvE22_KS4_eDB');
    $name = 'phormidolide_phmH_2_KS10_aMebOH' if($name eq 'phormidolide_phmH_2_KS10_0aOHaMebDOH');
    $name = 'phormidolide_phmF_5_KS8_bMeeDB' if($name eq 'phormidolide_phmF_5_KS8_eDBbMe');
    $name = 'phormidolide_phmI_1_KS13_bMeeDB' if($name eq 'phormidolide_phmI_1_KS13_eDBbMe');
    $name = 'phormidolide_phmF_3_KS6_adiMebOH' if($name eq 'phormidolide_phmF_3_KS6_adiMebDOH');
    $name = 'phormidolide_phmH_3_KS11_bOH' if($name eq 'phormidolide_phmH_3_KS11_bLOH');
    $name = 'phormidolide_phmI_5_KS_17_0bMeeDB' if($name eq 'phormidolide_phmI_5_KS17_0eDBbMe');
    $name = 'phormidolide_phmE_1_KS1_phosphoglycerate' if($name eq 'phormidolide_phmE_1_KS1_unusualStarter');
    $name = 'phormidolide_phmF_2_KS5_bOH' if($name eq 'phormidolide_phmF_2_KS5_bDOH');
    $name = 'phormidolide_phmF_4_KS7_bOH' if($name eq 'phormidolide_phmF_4_KS7_bDOH');
    $name = 'phormidolide_phmF_1_KS4_aMebOH' if($name eq 'phormidolide_phmF_3_KS4_aMebOH');
    $name = 'bryostatin_BryA_KS33_bDOH' if($name eq 'bryostatin_BryA_KS3_3bDOH');
    $name = 'chlorotonil_CtoD4_aMeshDB' if($name eq 'chlorotonil_CtoD4aMeshDB');
    $name = 'rhizoxins_RhiD1_K9_aMebDOH' if($name eq 'rhizoxins_RhiD1K9aMebDOH');
    $name = 'basiliskamides_P615_BasE1_KS2_aLMe_red' if($name eq 'basiliskamides_P615_BasE1_KS2_aLMered');
    $name = 'luminaolid_LumC2_KS8_aDMebDOH_bDOH' if($name eq 'luminaolid_LumC2_KS8_aDMebDOHbDOH');
    $name = 'luminaolid_LumC1_KS7_aLMebLOMe_aLMebLOH' if($name eq 'luminaolid_LumC1_KS7_aLMebLOMeaLMebLOH');
    $name = 'dorrigocin_migrastatin_MgsF_4KS7_aDMebLOH' if($name eq 'dorrigocin_migrastatin_MgsF4_KS7_aDMebLOH');
    $name = '9-methylstreptimidone_SmdI_KS55_aDMebketo' if($name eq '9-methylstreptimidone_SmdI_KS5_5aDMebketo');
    $name = 'misakinolide_MisC5_KS5_LOH_or_DPY' if($name eq 'misakinolide_MisC5_KS5_LOHorDPY');
    $name = 'etnangien_EtnH1_KS15_bimod_zDB' if($name eq 'etnangien_EtnH1_KS15_bimodzDB');
    $name = 'calyculin_CalE1_KS10_KS0ox' if($name eq 'calyculin_CalE_1_KS10_KS0ox');
    $name = 'tartrolon_TrtD_KS22_red' if($name eq 'tartrolon_TrtD_KS2_2red');
    $name = 'phormidolide_phmI_3_KS15_bOH' if($name eq 'phormidolide_phmI_3_KS15_bDOH');
        
    my $clade = $ln[$c];
    my $desc = join(' ', @ln[$c+1..$#ln]);
    #print "$name\n" unless(exists $leaf{$name});
    print join("\t", $name, $clade, $desc)."\n";
}
close $afh;
