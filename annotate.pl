#!/usr/bin/perl
use strict; use warnings;
use HTTP::Request;
use LWP::UserAgent;
use JSON;
use Getopt::Long;
my $help=''; my $outputfile=''; my $inputfile='';
#MANDATORY INPUT (-i) and OUTPUTFILE (-o). Will exit the program if not specified. 
GetOptions( 
	'o=s' => \$outputfile,
	'i=s' => \$inputfile,
	'h' => \$help) or die "Invalid args. Please specify input and output files.";
open (my $fh, "<", $inputfile);
open (my $fh_output, ">", $outputfile);
#display help message if -h invoked
if ($help){
	display_help();
	exit;
}
#if input and/or output file not specified, display help message and exit as they are mandatory
elsif ((!$outputfile) || (!$inputfile)){
	print "Invalid args. Please specify input and output filename";
	display_help(); 
	exit;
}
else{	
	#define tsv header
	print $fh_output "Chr\tPosition\tDepth\tNumber Reads\tPercentage Reads\tType Variant\tAllele Frequency\tExAC Gene(s)\t";
	print $fh_output "European (Non-Finnish)\tEuropean (Finnish)\tEast Asian\tSouth Asian\tLatino\tAfrican\tOther\n";
	my $i = 0; my @queries = (); 
	my @array=(); 
	while (my $row = <$fh>){
		#Ignore metadata
		if (!($row =~ /^#/) ){
			chomp $row;
			#perform bulk request in increments of 20 variants (limiting length of url)
			if ($i > 20){
				#pass the rows and formatted strings for json array to intiial subroutine 
				bulkify(\@array, \@queries);
				#reset increment back to zero and begin next set
				$i=0; @array=(); @queries=(); 
			}
			#gather the 4 variables necessary to append to the json array for ExAC bulk/variant/variant
			my ( $chr, $pos, $orig, $var ) = query_format($row);
			$i++;
			push @queries, "$chr-$pos-$orig-$var";
			push @array, $row;
		}
	}
	bulkify(\@array, \@queries);
}

close $fh; close $fh_output;

#run the initial subroutine that contains the necessary steps to produce a POST request
sub bulkify {
	my ( $array, $queries ) = @_;
	my $json_string = encode_json($queries);
	my ( $decoded ) = call($json_string);
	my ( @formatted_rows )= bulk($array, $decoded, $queries);
}

#output the correct format required for each entry of the json array
sub query_format{
	my ( $row ) = @_;
	my @splitrow = split("\t", $row);
	my ( $most_deleterious ) = $splitrow[4] =~ /^(.+?),?/;
	return($splitrow[0], $splitrow[1], $splitrow[3], $most_deleterious)
}

#generate the percentage from the fraction
sub percent_variation {
	my ( $proportion_var )= @_;
	my $percent_var = sprintf("%.4f", $proportion_var)*100; 
	return $percent_var;
}
sub bulk {
	my ( $ba, $decoded, $quer ) = @_;
	my @bulk_array = @{$ba};
	my @queries = @{$quer};
	my ( $freq, $pop_homs, $genes );
	for (my $i=0; $i < scalar(@bulk_array); $i++){
		#split each row into individual columns
		my @splitrow = split("\t", $bulk_array[$i]);
		#Capture the type of variant via regex. Most deleterious type listed first
		my ( $type ) = $splitrow[7] =~ /TYPE=([^,]*)/;
		#Capture read Depth. Most deleterious listed first. 
		my ( $depth ) = $splitrow[7] =~ /DP=([^,;]*)/;
		#Capture variant count. Most deleterious listed first. 
		my ( $sup_variant ) = $splitrow[7] =~ /AC=([^,;]*)/;
		#Capture proportion and percentage of allele frequency. Most deleterious listed first.
		my ( $proportion_var ) = $splitrow[7] =~ /AF=([^,;]*)/;
		my $percent_var = percent_variation($proportion_var);
		if (exists $decoded->{$queries[$i]}->{allele_freq}){
	 		( $freq, $pop_homs, $genes )= additional_information($decoded, $queries[$i]);	
	 	}
	 	else{
	 		$freq=$genes=".";
	 		$pop_homs = ".\t" x 7;
	 	}
	 	#print everything in format according to tsv header
		print $fh_output $splitrow[0]."\t".$splitrow[1]."\t".$depth."\t";
		print $fh_output $sup_variant."\t".$percent_var."\t".$type."\t".$freq."\t".$genes."\t".$pop_homs."\n";
		
	}

}

sub call{
	my ( $json_string ) = @_;
	my ( $freq, $pop_homs, $genes );
	#Retrieve Allele freq of variant from ExAC API (Bulk query)
	my $str = "http://exac.hms.harvard.edu/rest/bulk/variant/variant";
	my $Useragent = LWP::UserAgent->new;
	my $json = $json_string;
	#bulk requires post request with json array as body
	my $req = HTTP::Request->new('POST', $str);
	#specify header
	$req -> header('Content-Type'=>'application/json');
	#give the json array (as string) as the body
	$req->content( $json );
	my $returned = $Useragent->request( $req );
	#decode the newly obtained hash from POST and select only the returned content portion
	my $decoded = decode_json($returned->{_content});
	return $decoded;
}
sub additional_information{
	my ( $decoded, $index ) = @_;
	my @array_genes=();
	my @pop_homs=();
	my ( $freq, $pop_homs, $genes );
	#give frequency of allele
	if (exists $decoded->{$index}->{allele_freq}){
		$freq = sprintf("%.6f", $decoded->{$index}->{allele_freq});
	}
	#should it exist, update the genes associated with variant from ExAC
	if (exists $decoded->{$index}->{genes}){
		my @array_genes = @{$decoded->{$index}->{genes}};
		if (scalar(@array_genes) != 0){
			foreach (@array_genes){
				
				$genes .= $_.",";
			}
		}
		else{
			$genes = ".";
		}
		$genes =~ s/,$//;
	}
	#gather the individual population numbers 
	if (exists $decoded->{$index}->{pop_homs}){
		$pop_homs .= $decoded->{$index}->{pop_homs}->{'European (Non-Finnish)'}."\t";
		$pop_homs .= $decoded->{$index}->{pop_homs}->{'European (Finnish)'}."\t";
		$pop_homs .= $decoded->{$index}->{pop_homs}->{"East Asian"}."\t";
		$pop_homs .= $decoded->{$index}->{pop_homs}->{"South Asian"}."\t";
		$pop_homs .= $decoded->{$index}->{pop_homs}->{"Latino"}."\t";
		$pop_homs .= $decoded->{$index}->{pop_homs}->{"African"}."\t";
		$pop_homs .= $decoded->{$index}->{pop_homs}->{"Other"};
		
	}
	else {
		print "doesnt exists\n";
		$pop_homs = ".\t.\t.\t." x 7;
		$pop_homs =~ s/,$//;
	}
	return ($freq, $pop_homs, $genes);
}
#display help message with -h argument
sub display_help {
	print "Usage: perl annotate.pl [OPTION] [filename]\n";
	print "\t-i, input VCF file\n\t-o, output tsv file\n";
	print "This is a short perl script for annotating a VCF and extracting additional\n";
	print "additional information from ExAC\n";

}