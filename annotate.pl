#!/usr/bin/perl
use strict; use warnings;
use HTTP::Request;
use LWP::UserAgent;
use JSON;
use Getopt::Long;
my $help=''; my $outputfile=''; my $inputfile='';

GetOptions( 
	'o=s' => \$outputfile,
	'i=s' => \$inputfile,
	'h' => \$help) or die "Invalid args. Please specify input and output files.";
open (my $fh, "<", $inputfile);
open (my $fh_output, ">", $outputfile);
if ($help){
	display_help();
	exit;
}
elsif ((!$outputfile) || (!$inputfile)){
	print "Invalid args. Please specify input and output filename";
	display_help(); 
	exit;
}
else{	
	print $fh_output "Chr\tPosition\tDepth\tNumber Reads\tPercentage Reads\tType Variant\tAllele Frequency\tExAC Gene(s)\t";
	print $fh_output "European (Non-Finnish)\tEuropean (Finnish)\tEast Asian\tSouth Asian\tLatino\tAfrican\tOther\n";
	my $i = 0; my @queries = (); 
	my @array=(); 
	while (my $row = <$fh>){
		#Ignore metadata
		if (!($row =~ /^#/) ){
			chomp $row;
			if ($i > 20){
				bulkify(\@array, \@queries);
				$i=0; @array=(); @queries=(); 
			}
			my ( $chr, $pos, $orig, $var ) = query_format($row);
			$i++;
			push @queries, "$chr-$pos-$orig-$var";
			push @array, $row;
		}
	}
	bulkify(\@array, \@queries);
}

close $fh; close $fh_output;

sub bulkify {
	my ( $array, $queries ) = @_;
	my $json_string = encode_json($queries);
	my ( $decoded ) = call($json_string);
	my ( @formatted_rows )= bulk($array, $decoded, $queries);
}


sub query_format{
	my ( $row ) = @_;
	my @splitrow = split("\t", $row);
	my ( $most_deleterious ) = $splitrow[4] =~ /^(.+?),?/;
	return($splitrow[0], $splitrow[1], $splitrow[3], $most_deleterious)
}


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
	my $req = HTTP::Request->new('POST', $str);
	$req -> header('Content-Type'=>'application/json');
	$req->content( $json );
	my $returned = $Useragent->request( $req );
	my $decoded = decode_json($returned->{_content});
	return $decoded;
}
sub additional_information{
	my ( $decoded, $index ) = @_;
	my @array_genes=();
	my @pop_homs=();
	my ( $freq, $pop_homs, $genes );
	if (exists $decoded->{$index}->{allele_freq}){
		$freq = sprintf("%.6f", $decoded->{$index}->{allele_freq});
	}
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
sub display_help {
	print "Usage: perl annotate.pl [OPTION] [filename]\n";
	print "\t-i, input VCF file\n\t-o, output tsv file\n";
	print "This is a short perl script for annotating a VCF and extracting additional\n";
	print "additional information from ExAC\n";

}