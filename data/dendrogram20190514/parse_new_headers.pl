#!/bin/env perl

use strict;
use warnings;

while(<>){
    chomp;
    if($_=~ m/^>(.+)/){
	my ($genome, $k, $gloc, $strand, $locustag, $ann, $prot, $geneloc, $prot_ks, $kks_ploc) = split(/\|/, $1);
	my ($k_ks, $ploc) = split(/\//, $kks_ploc);

	if($k_ks eq 'KS55'){
	    $k_ks = 'KS5';
	}
	if($genome eq 'HE804045' && $k eq 'c19'){
	    if($locustag eq 'BN6_45860' && $prot_ks eq 'KS1'){
		$k_ks = 'KS1';
	    }elsif($locustag eq 'BN6_45870' && $prot_ks eq 'KS1'){
		$k_ks = 'KS2';
	    }elsif($locustag eq 'BN6_46060' && $prot_ks eq 'KS1'){
		$k_ks = 'KS3';
	    }elsif($locustag eq 'BN6_46170' && $prot_ks eq 'KS1'){
		$k_ks = 'KS4';
	    }elsif($locustag eq 'BN6_46170' && $prot_ks eq 'KS2'){
		$k_ks = 'KS5';
	    }elsif($locustag eq 'BN6_46170' && $prot_ks eq 'KS3'){
		$k_ks = 'KS6';
	    }elsif($locustag eq 'BN6_46170' && $prot_ks eq 'KS4'){
		$k_ks = 'KS7';
	    }elsif($locustag eq 'BN6_46180' && $prot_ks eq 'KS1'){
		$k_ks = 'KS8';
	    }elsif($locustag eq 'BN6_46180' && $prot_ks eq 'KS2'){
		$k_ks = 'KS9';
	    }elsif($locustag eq 'BN6_46180' && $prot_ks eq 'KS3'){
		$k_ks = 'KS10';
	    }
	}

	if($genome eq 'CP002399' && $k eq 'c9'){
	    if($locustag eq 'ML5_4521' && $prot_ks eq 'KS1'){
		$k_ks = 'KS1';
	    }elsif($locustag eq 'ML5_4521' && $prot_ks eq 'KS2'){
		$k_ks = 'KS2';
	    }elsif($locustag eq 'ML5_4522' && $prot_ks eq 'KS1'){
		$k_ks = 'KS3';
	    }elsif($locustag eq 'ML5_4522' && $prot_ks eq 'KS2'){
		$k_ks = 'KS4';
	    }elsif($locustag eq 'ML5_4522' && $prot_ks eq 'KS3'){
		$k_ks = 'KS5';
	    }elsif($locustag eq 'ML5_4522' && $prot_ks eq 'KS4'){
		$k_ks = 'KS6';
	    }elsif($locustag eq 'ML5_4534' && $prot_ks eq 'KS1'){
		$k_ks = 'KS7';
	    }
	}

	if($genome eq 'CP002162' && $k eq 'c12'){
	    if($locustag eq 'Micau_3883' && $prot_ks eq 'KS1'){
		$k_ks = 'KS1';
	    }elsif($locustag eq 'Micau_3895' && $prot_ks eq 'KS1'){
		$k_ks = 'KS2';
	    }elsif($locustag eq 'Micau_3895' && $prot_ks eq 'KS2'){
		$k_ks = 'KS3';
	    }elsif($locustag eq 'Micau_3895' && $prot_ks eq 'KS3'){
		$k_ks = 'KS4';
	    }elsif($locustag eq 'Micau_3895' && $prot_ks eq 'KS4'){
		$k_ks = 'KS5';
	    }elsif($locustag eq 'Micau_3896' && $prot_ks eq 'KS1'){
		$k_ks = 'KS6';
	    }elsif($locustag eq 'Micau_3896' && $prot_ks eq 'KS2'){
		$k_ks = 'KS7';
	    }
	}

	if($genome eq 'misakinolide'){
	    if($k_ks eq 'KS8'){
		$k_ks = 'KS7';
	    }elsif($k_ks eq 'KS9'){
		if($locustag eq 'MisD2'){
		    $k_ks = 'KS8';
		}
	    }
	}

	if($genome eq 'basiliskamides'){
	    if($kks_ploc =~ m/^KS3\/36/){
		$k_ks = 'KS4';
	    }
	}

	if($genome eq 'disorazole'){
	    if($ann eq 'DszC1' && $k_ks eq 'KS2'){
		$k_ks = 'KS10';
	    }
	}

	if($genome eq 'etnangien'){
	    if($ann eq 'EtnE3' && $k_ks eq 'KS3'){
		$k_ks = 'KS7';
	    }
	}

	if($genome eq 'misakinolide'){
	    if($ann eq 'MisD2'){
		if($k_ks eq 'KS8'){
		    $k_ks = 'KS7';
		}elsif($k_ks eq 'KS9'){
		    $k_ks = 'KS8';
		}
	    }
	}

	if($genome eq 'myxovirescin'){
	    if($ann eq 'TaP1'){
		$k_ks = 'KS12';
	    }
	}

	if($genome eq 'nosperin'){
	    if($ann eq 'NspC5'){
		$k_ks = 'KS8';
	    }
	}

	if($genome eq 'rhizoxins'){
	    if($ann eq 'RhiE1'){
		$k_ks = 'KS12';
	    }
	}
	
	if($genome eq 'bongkrekic_acid'){
	    if($kks_ploc =~ m/^KS11\/33/){
		$k_ks = 'KS1';
	    }
	}

	if($genome eq 'calyculin'){
	    if($ann eq 'CalE6'){
		$k_ks = 'KS15';
	    }elsif($ann eq 'CalE7'){
		$k_ks = 'KS16';
	    }elsif($ann eq 'CalF1'){
		$k_ks = 'KS17';
	    }elsif($ann eq 'CalF2'){
		$k_ks = 'KS18';
	    }elsif($ann eq 'CalF3'){
		$k_ks = 'KS19';
	    }elsif($ann eq 'CalF4'){
		$k_ks = 'KS20';
	    }elsif($ann eq 'CalF5'){
		$k_ks = 'KS21';
	    }elsif($ann eq 'CalG1'){
		$k_ks = 'KS22';
	    }elsif($ann eq 'CalG2'){
		$k_ks = 'KS23';
	    }elsif($ann eq 'CalG3'){
		$k_ks = 'KS24';
	    }elsif($ann eq 'CalG4'){
		$k_ks = 'KS25';
	    }elsif($ann eq 'CalH1'){
		$k_ks = 'KS26';
	    }elsif($ann eq 'CalH2'){
		$k_ks = 'KS27';
	    }elsif($ann eq 'CalH3'){
		$k_ks = 'KS28';
	    }elsif($ann eq 'CalI1'){
		$k_ks = 'KS29';
	    }
	}

	if($genome =~ m/^2hg4_chainA_EryKS5/){
	    $genome = 'OUTGROUP_2hg4';
	    $k = 'chainA_EryKS5';
	    $k_ks = 'KS1';
	}elsif($genome =~ m/^2qo3_chainA_EryKS3/){
	    $genome = 'OUTGROUP_2qo3';
	    $k = 'chainA_EryKS3';
	    $k_ks = 'KS1';
	}
	
	if($k eq 'acc'){
	    print '>'.join("|", $genome, $k_ks)."\n";
	}else{
	    print '>'.join("|", join('|', $genome, $k), $k_ks)."\n";
	}
    }else{
	print "$_\n";
    }
}
