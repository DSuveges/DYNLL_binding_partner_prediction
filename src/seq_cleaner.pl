#!/usr/bin/env perl
#version 3.1 last edited: 10.08.2010

#This program reads a multifasta file
#perform a disorder predictioon -> results stored in a binary number, where 1 
#corresponds to a disordered residue and 0 corresponds to a ordered residue

#The program also perform a check for transmembrane or extracellular domains of proteins.
#these domains (as they are excluded from the search) masked with an X

# Requirements:
    # *Internet connection
    # *LPW perl library
    # *iupred in PATH
    
# Usage:
# $./seq_cleaner.pl -i <inputfile.fasta> -o <outputfile>

# !!Please notice that the program is fit to the header of the normal UniProt fasta file.
# In case of other formats please modify the regexp at line "28"

use Getopt::Std;
getopt('i:o:');

open (INFILE, "<", $opt_i);
open (OUTFILE, ">", $opt_o);
foreach ( <INFILE> ){
	if ($_ =~ /^>/){

        # Here you can fit the program to your own fasta header type
        $_ =~ /^>UniRef90_(.6)+\s(.+?)n=/;
		print "Processing - $1\t$2\n";

		if ($sequence) {
			$sequence  =~ s/\s//g;
			$cytosolic =  &cytosolic($protein_id,$protein_name,$sequence);
			$disorder  =  &disorder($sequence);
			print OUTFILE "$protein_id\t$protein_name\t$cytosolic\t$disorder\n";
		}

        $sequence       = '';
		$protein_id     = $1;
		$protein_name   = $2;
	}
	else {
		$sequence = $sequence.$_;
	}    
}

# Necessary for the last sequence
if ($sequence) {
    $sequence  =~ s/\s//g;
    $sequence  =~ s/\s//g;
	$cytosolic =  &cytosolic($protein_id,$protein_name,$sequence);
	$disorder  =  &disorder($sequence);
	print OUTFILE "$protein_id\t$protein_name\t$cytosolic\t$disorder\n";
	($sequence, $protein_id,$protein_name) = '';
}
close INFILE;
close OUTFILE;

# This subroutine masks transmembrane and extracellular domains
sub cytosolic {

    my $Prot_ID  = shift(@_);
	my $name     = shift(@_);
    my $sequence = shift(@_);
	my $new_name = '';
    my $newseq   = '';
    my $segment  = '';
    my $length   = '';
    my $tag      = '0';
	
    use LWP::UserAgent;
	use LWP::Simple;
	use URI;
	
    my $url = 'http://www.uniprot.org/uniprot/'.$Prot_ID.'.txt';
	my $entry = get $url;
	# my $cytosol_tag = '';
	my $tm_tag = '0';
	my @regions = '';
	
    my @entry = split(/\n/,$entry);

	foreach my $line (@entry){
		if ($line =~ /^FT\s+TRANSMEM\s+(\d+)\s+(\d+)\s+/){
            push (@regions, "$1-$2");
            $new_name = $name." Transmembrane protein";

        }
        if ($line =~ /^FT\s+TOPO_DOM\s+(\d+)\s+(\d+)\s+Extracellular/){
            push (@regions, "$1-$2");
		    $new_name = $name." Transmembrane protein";

        }
	}
	
	my $newsequence = '';
	my($lower,$upper, $old_upper) = '';

	
	my $upper_old   = "0";
	my $lower_old   = length($sequence);
	
	shift @regions;
	foreach $regio (@regions){
		($lower,$upper)    = split(/-/,$regio);
		# print "$lower ... $upper\n";
		my $segment = substr($sequence,$upper_old,$lower-$upper_old);
		$newsequence = $newsequence.$segment;
		for(1..$upper-$lower){$newsequence .= "X"}
		$upper_old = $upper;
	}
	if ($upper <= length $sequence){$newsequence .= substr($sequence,$upper)}
	
	return $newsequence;

}

# This subroutine performs the disorder prediction.
sub disorder {
    #print a dot for each calling, that is equal the number of processed sequences!
    print ".";

    my $seq    = shift(@_);
	#A temporary fasta sequence file is required
   	open (TEMP, ">", "temp.fasta");
   	print TEMP ">header\n$seq";
   	close TEMP;
   	#run iupred for long disordered segment
   	my @iupred_prediction = qx(iupred temp.fasta long); 
   
    my $disorder = '';
    
	#the iupred prediction is transfered into a binary number, 
	#where 0 corresponds to ordered, 1 to disordered residue.
    #A residue is considered to be disordered if the disorder tendency >= 0.5
    foreach my $line (@iupred_prediction){
		unless ($line =~ /^#/){	#Ignoring header lines from the output
			$line   =~ /^\s+(\d+)\s+[A-Z,a-z]{1}\s+([0]\.[0-9]{4})/;
			if ($2 >= 0.5){
			    $disorder = $disorder."1";#disordered
			}
			else {
			    $disorder = $disorder."0";#ordered
			}
		}
	}

    return $disorder;
}

