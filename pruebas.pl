#!/usr/bin/perl -w

use BIO::Seq;
use BIO::DB::SoapEUtilities;

my $fac = Bio::DB::SoapEUtilities->new();
my $esearch = $fac->esearch(-db => 'nucleotide',
			    -term => 'H1V1 and CCR5 and Brazil');#'Human papillomavirus type 2 AND Complete Genome',
my $results = $esearch->run(-auto_adapt => 1, -rettype => 'fasta');
#print "Query translation: ", $fac->get_query_translation(), "\n";
#print "Count = ", $results->count(), "\n";
#my @ids = $fac->get_ids();
#foreach $i (@ids){
#   print $i,"\n";
#}
