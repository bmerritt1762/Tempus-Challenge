# Tempus-Challenge
Short Script to Analyze VCF and Annotate. Additional information for each variant is gathered via ExAC. 
Variants are queried in bulk with ExAC API to speed up computation time. 
Only the most deleterious possibility is selected from each variant entry.

Language: Perl 5, version 22
Modules used: JSON, LWP, HTTP, Getopt

Usage: perl annotate.pl [OPTION] [filename]
	-i, input VCF file
	-o, output tsv file

Output header:

1. Chromosome
2. Position of variant
3. Depth of sequence coverage
4. Number of reads supporting variant
5. % of Reads supporting variant versus reference reads
6. Type of variant (del, snp, ins, etc)
7. Allele Frequency of variant in population
8. ExAC listed genes associated with variant
9-15. Population counts (pop_homs) 


