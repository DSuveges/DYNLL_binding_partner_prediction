#!/usr/bin/env perl
#version 3.0 last edited: 22.02.2011

# This program searches potential DYNLL binding sites in the disordered regions of proteins
# Uses the output of the seq_cleaner.pl (v3.1)
#Exclude the fully ordered segemnts, segments that contain transmembrane or extracellular residues, and segments that  do not contain gutamine at 0th position.

# usage:
# $./predictor.pl -s <scorefile> -i <infile> -o <outfile> -d <score_distrib>

    # <scorefile> - the output of the motif2score.pl program
    # <infile>    - the output of the seq_cleaner.pl program
    # <outfile>   - the detailed output, with the binding sequences, scores, positions etc...
    # <score_distrib> - distribution of scores.

use Getopt::Std;
getopt('s:i:o:d:');

our %score_hash = '';

#this variable should be modified in scrambled score matrices
our $core_position = '5';
our %score_distribution = '';

for ($j = "100"; $j <= "480"; $j++){
    $score_distribution{$j} = "0";
}

#score reading 
open (SCOREFILE, "<", "$opt_s");#open score file
foreach ( <SCOREFILE> ){
    $_ =~ s/\n//;
   	my ($position, $aminoacid, $score) = split (/\t/,$_);
    $score_hash{$position}{$aminoacid} = $score;
}
close SCOREFILE;

open (INFILE, "<", $opt_i);
open (OUTFILE, ">", $opt_o);
foreach (<INFILE>){
    $_ =~ s/$\n//;
    my ($prot_ID,$name,$prot_seq,$disorder) = split(/\t/,$_);
    &scoring($prot_seq,$disorder,$name,$prot_ID);
}
close INFILE;
close OUTFILE;

open (DISTRIB, ">", $opt_d);
foreach $scores (keys %score_distribution){
    print DISTRIB "$scores\t$score_distribution{$scores}\n";
}
close DISTRIB;


sub scoring {
    my $seq      = $_[0];
    my $disord   = $_[1];
    my $name     = $_[2];
	my $ID		 = $_[3];
    my $distance = '';
    
    #creating all possible eight residue long segment of the sequence.
	for ( my $position = "0"; $position + 8 <= length($seq); $position++){

		my $segment = substr ($seq, $position, 8);
        my $score = '';
        
        #considering sequences where there is glutamine in the appropriate position
        if ((substr($segment,$core_position,1) eq "Q") && ($segment !~ /x/i)){

			$score = &scoring_routine($segment);
            my $disord_number  = substr ($disord, $position, 8);
            # print "$segment\t$disord_number\t$score\n";            
            #considering onli disordered seuences
			if ($disord_number > "0"){
	            if ($score > "200" ) {
                    print OUTFILE "$ID\t$name\t$position\t$segment\t$score\n";
                }
			    $score_distribution{$score}++;
			}
		}   
	}
}

sub scoring_routine {
    my $binding_peptide = shift( @_ );
    #print "$binding_peptide\n";
    my $aminoacid = '';
    my $position = 0;
    my $score = 0;
    
    #We check the score for each amino acid within the global 
	#"%score_hash" -> the scores are summarized for each position!
    for ($position = 0; $position <= 7; $position++){
        $aminoacid = substr ($binding_peptide,$position,1);
        $score = $score + $score_hash{$position}{$aminoacid};
    }
    
    return $score;
}
