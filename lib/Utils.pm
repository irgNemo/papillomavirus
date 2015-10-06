package Utils;

use LWP::Simple;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw (getNCBIIDs);

sub getCompleteGenomeNCBIIDs {

my ($database, $organism) = (@_);
$organism="%22".$organism."%22";
$organism
my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=".$database."&term=


#%22Human%20papillomavirus%20type%202%22[Organism]%20OR%20%22Human%20Papillomavirus%20type%202%22[All%20Fields]%29%20AND%20%22complete%20genome%22[All%20Fields]";
my $content = get ($url);


}


my @matches = map {s/\<\/Id\>//g; $_;} map {s/\<Id\>//g; $_; } grep { /\<Id\>[0-9]+\<\/Id\>/ } split(/\n/,$content);

#for(my $i=0;$i<=$#matches;$i++){
#	$url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=".$matches[$i]."&rettype=gb";
#	$filename=$matches[$i].".gb";
#	getstore ($url,$filename);
#	#print $url;
#} 
