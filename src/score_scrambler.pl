#!/usr/bin/perl
#version 1.0 last edited: 10.01.2010

# This program reads the scorefile created by motif2score.pl
# Randomize the positions of the scores 1060 times. 
# (Required to get ~1000 different scorefiles)
# Writes the output scorefiles named after the consensus sequence.

# For randomization, the program uses the Fisher Yates Shuffle algorithm:
# http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle

# usage:
# $./score_scrambler.pl -s <scorefile>


use Getopt::Std;
getopt('s:');

@array = (0..7);

our %score_hash = '';

open (SCOREFILE, "<", "$opt_s");
foreach ( <SCOREFILE> ){
    $_ =~ s/\n//;
   	my ($position, $aminoacid, $score) = split (/\t/,$_);
 	$score_hash{$position}{$aminoacid} = $score;
}
close SCOREFILE;

#The number of randomization can be set here:
for (1...60){
    my $back = &fisher_yates_shuffle(\@array);
    my $peptide = &peptide($back);
    &output($peptide,$back);    
}

sub fisher_yates_shuffle {
    my $array = shift(@_);
    my $i;
    
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
    return \@array;
}

sub peptide {
    my $array = shift(@_);
    my $consensus = "VSRGTQTE";
    my $new_sequence = '';
    
    for ($i = "0"; $i <= "7" ; $i++){
        my $aa = substr($consensus,${$array}[$i],1);
        $new_sequence = $new_sequence.$aa;
    }
    
    return $new_sequence;
}

sub output {
    my $peptide = shift(@_);
    my $array_ref = shift(@_);
    open (OUTPUT,">$peptide");
    for ($i = "0"; $i <= "7"; $i++) {
        foreach $aminoacid (sort keys %{$score_hash{$array_ref->[$i]}}){
            print OUTPUT "$i\t$aminoacid\t$score_hash{$array_ref->[$i]}{$aminoacid}\n";
        }
    }
    close OUTPUT
}