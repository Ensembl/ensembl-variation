#!/usr/bin/perl

use strict;
use warnings;

#### Check genomic, coding & non-coding HGVS strings give variation features which return the same HGVS strings
#### This is gene-annotation dependant


use Test::More;

use Data::Dumper;

use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

BEGIN {
    use_ok('Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor');
}
my $DEBUG = 0;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_all("$Bin/test.ensembl.registry");

my $variationfeature_adaptor    = $reg->get_adaptor('human', 'variation', 'variationfeature');
my $transcript_adaptor          = $reg->get_adaptor('human', 'core',      'transcript');
my $transcript_variation_adaptor= $reg->get_adaptor('human', 'variation', 'transcriptvariation');


## TEST DATA:  hgvs_genomic,  variant_allele, hgvs_[non]coding,  variant_allele, hgvs_protein,  test_description 

my @test_input = (    ["2:g.46746465G>A",    
		       "A", 
		       "ENST00000522587.1:c.-101-6514C>T", 
		       "T",
		       "",           
		       "substitution, coding intron - downstream"
		      ], 

		      ["12:g.102061070_102061071insT",
		       "T",
		       "ENST00000360610.2:c.2336-440_2336-439insT", 
		       "T",
		       "",
		       "insertion, coding intron downstream"
		      ],
		      ["4:g.41993003G>A",    
		       "A", 
		       "ENST00000264451.6:c.109+226G>A",   
		       "A",
		       "",           
		       "substitution, coding intron upstream"], 
		      ["2:g.46739156C>T",    
		       "T",
		       "ENST00000306448.4:c.*14G>A",       
		       "A",
		       "",           
		       "substitution 3' UTR" ], 
		      ["2:g.46746256G>A",    
		       "A", 
		       "ENST00000306448.4:c.-274C>T",      
		       "T",
		       "",           
		       "substitution, 5' UTR" 
		      ],  
		      ["2:g.46746507delCinsACAA",
		       "ACAA",  
		       "ENST00000524249.1:n.775+16445delGinsTTGT",   
		       "TTGT",
		       "", 
		       "delins, non-coding "
		      ],  
		      ["2:g.46739212C>G",    
		       "G", 
		       "ENST00000522587.1:c.639G>C",  
		       "C", 
		       "ENSP00000428141.1:p.Met213Ile",      
		       "substitution, non_syn"
		      ],

		      ["12:g.102056227G>A",
		       "A",
		       "ENST00000360610.2:c.2049G>A",
		       "A",
		       "ENST00000360610.2:c.2049G>A(p.=)",
		       "substitution synonymous"
		      ],
		      ["2:g.46739488G>A",    
		       "A", 
		       "ENST00000306448.4:c.363C>T",  
		       "T", 
		       "ENST00000306448.4:c.363C>T(p.=)",    
		       "substitution, synonymous"
		      ],     
		      ["7:g.143557504delT",  
		       "-", 
		       "ENST00000355951.2:c.1964delA", 
		       "-",
		       "ENSP00000348220.2:p.Gln655ArgfsX17",  
		       "deletion, fs"
		      ], 		      
		      ["7:g.7680048A>G",
		       "G",
		       "ENST00000223129.4:c.2T>C",
		       "C",
		       "ENSP00000223129.4:p.Met1?",
		       "substitution,  start loss"
		      ],
		      ["X:g.131215393dupA",
		       "A",
		       "ENST00000298542.4:c.905+997dupT", 
		       "T",
		       "",
		       "duplication, intronic rc transcript"
		      ], 
		      ["6:g.6649978_6649980dupAGG",
		       "AGG",
		      "ENST00000230568.3:c.405+68_405+70dupAGG",
		      "AGG",
		       "",
		       "duplication, intronic - long"
		      ],
		      ["22:g.20920895_20920939dupCCACAGCCTCCGCCCTCCCAGGCTCTGCCCCAGCAGCTGCAGCAG",
		       "CCACAGCCTCCGCCCTCCCAGGCTCTGCCCCAGCAGCTGCAGCAG",
		       "ENST00000292733.7:c.832_876dupCCACAGCCTCCGCCCTCCCAGGCTCTGCCCCAGCAGCTGCAGCAG",
		       "CCACAGCCTCCGCCCTCCCAGGCTCTGCCCCAGCAGCTGCAGCAG",
		       "ENSP00000292733.7:p.Pro278_Gln292dup",		       
		       "insertion, peptide duplication"
		      ],
		      ["6:g.31997361delC",
		       "-",
		       "ENST00000435363.2:c.3695delC",           
		       "-",
		       "ENSP00000415941.2:p.Ser1232X",
		       "deletion, stop gained"
		      ],		      
		      ["4:g.130032945A>G",
		       "G",
		       "ENST00000281146.4:c.599A>G",
		       "G",
		       "ENSP00000281146.4:p.X200TrpextX2",
		       "substitution, stop lost"
		      ],		    
		      ["1:g.154140413_154140415delTTA",
		       "-", 
		       "ENST00000368530.2:c.856_858delTAA",
		       "-",
		       "ENSP00000357516.2:p.X286delextX56",
		        "deletion, stop loss"
		      ],
		      [ "19:g.48836478_48836480delAGG", 
			"-",
			"ENST00000293261.2:c.1376_1378delCCT", 
			"-",
			"ENSP00000293261.2:p.Ser459del", 
			"deletion, inframe codon loss"
		      ],
		      ["2:g.46731836A>G",    
		       "G",  
		       "",                          
		       "",
		       "",
		       "substitution, downstream" 
		      ], 
		      ["2:g.46747460C>G",    
		       "G", 
		       "",
		       "",
		       "",
		       "substitution, upstream" 
		      ], 
		      ["2:g.46732522G>A",    
		       "A", 
		       "ENST00000524249.1:n.776-14552C>T", 
		       "T",
		       "",           
		       "substitution, noncoding intron"
		      ], 
);

foreach my $line(@test_input){

    if($DEBUG==1){     print "\n\n\nStarting $line->[0]/ $line->[2]\n";}

    ## create variation feature from hgvs string
    my $variation_feature_g = $variationfeature_adaptor->fetch_by_hgvs_notation($line->[0] );
    test_output($variation_feature_g, "genomic", $line, $line->[1] );

    if( $line->[2] =~/\w+/){
	my $variation_feature_c = $variationfeature_adaptor->fetch_by_hgvs_notation($line->[2] );    
	test_output($variation_feature_c, "transcr" , $line, $line->[3]);
    }
}

sub test_output{

    my $variation_feature = shift;
    my $input_type        = shift;
    my $line              = shift;
    my $allele            = shift;


    ##### genomic
    my $hgvs_genomic      = $variation_feature->get_all_hgvs_notations("", "g");

    ok( $line->[0] eq $hgvs_genomic->{$allele} , "$input_type -> genomic - $line->[5]  [$line->[0]] ");
    if($DEBUG==1){print "TEMP: $line->[0]\t$hgvs_genomic->{$line->[1]}\n";}					   

    
    ##### transcript level - transcript to be supplied as ref feature  - alt allele may be complimented wrt genomic reference  
    if( $line->[2] =~/\w+/){

	my $transcript_name = (split/\./, $line->[2])[0];
        my $transcript      = $transcript_adaptor->fetch_by_stable_id( $transcript_name);   

         my $hgvs_coding      = $variation_feature->get_all_hgvs_notations($transcript, "c");

	ok( $line->[2] eq $hgvs_coding->{$allele} , "$input_type -> transcr - $line->[5] [$line->[2]]");
         if($DEBUG==1){	print "TEMP: $line->[2]\t$hgvs_coding->{$line->[1]} [allele:$line->[1]]\n";}


        if( $line->[4] =~/\w+/){
	    ##### protein level - transcript to be supplied as ref feature - alt allele may be complimented
	    
	    #my $transcript_variation = $transcript_variation_adaptor->fetch_all_by_VariationFeatures( [$variation_feature], [$transcript]);
	    my $transcript_variations = $variation_feature->get_all_TranscriptVariations();
	    foreach my $transcript_variation (@{$transcript_variations}){

		if($transcript_variation->transcript->stable_id() eq $transcript_name){
		    
		    my $hgvs_protein      = $variation_feature->get_all_hgvs_notations($transcript, 
										       "p",
										       $transcript_name,
										       "",
										       $transcript_variation
			);
		    ok( $line->[4] eq $hgvs_protein->{$allele} , "$input_type -> protein - $line->[5] [$line->[4]]");
		    if($DEBUG==1){   print "TEMP: $line->[4]\t$hgvs_protein->{$line->[3]}\n";}
		}
	    }
	}
    }
}

done_testing();

