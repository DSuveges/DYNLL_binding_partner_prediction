#!/usr/bin/env perl
#Version.: 2.0
#Last modified: 10.11.2010

#This program reads a list of DYNLL binding sequences (-i), and 
#creates a three dimensional matrix with the probability of each aminoacid in each position.
#the matrix is written into the output file (-o).

#Usage:
#$./motif2score.pl -i <infile> -o <outfile>

use Getopt::Std;
getopt('i:o:');

my $aa = '';
my $seq_number = "0";

open (INFILE, "<", "$opt_i");
foreach ( <INFILE> ){
	$_ =~ s/\s//;
	my $pos = "0";
	
	while ($pos  <  length $_ ){
		$aa = substr($_, $pos, "1");
		$pos++;
		unless($hash{$pos} =~ /$aa/i){
			$hash{$pos} .= $aa;
		}
	}
    
	print "\n";
	$pos = "0";
	$seq_number ++;
}
close INFILE;

open (OUTFILE, ">", "$opt_o");
for $position ( sort keys %hash) {
	for $aminoacid (sort keys %{$hash{$position}}){
		$hash{$position}{$aminoacid} = int ( 100 * $hash{$position}{$aminoacid} / $seq_number );
		print OUTFILE "$position\t$aminoacid\t$hash{$position}{$aminoacid}\n";
    }
}
close OUTFILE;
