package Utils;
use strict;
use warnings;
use Exporter;
use LWP::Simple;
use Data::Dumper;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Utils);
our @EXPORT = qw(getNCBICompleteGenomebyID createRandomSets);


sub getNCBICompleteGenomebyID {
	my ($id,$path) = @_ or die "Wrong parameter number in getNCBICompleteGenomebyOrganism function";
	my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=".$id."&rettype=gb";
	my $filename=$path.$id.".gb";
	getstore ($url,$filename);
}

1;
