package Transformer;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Transformer);
our @EXPORT = qw(statisticsPerTag);

sub SeqIOToHash{
	my ($seqio_obj, $tagsToVerify, $primaryTag, $secondaryTag) = @_ or die "Wrong parameter number on Statistic::statisticPerTag function";
	my %sequences;
	while(my $seq = $seqio_obj->next_seq){
		my %tagsFound;
		foreach my $featureObj ($seq->get_SeqFeatures){
			foreach my $innerTag ($featureObj->get_all_tags){
				if (($secondaryTag eq $innerTag) and ($featureObj->primary_tag eq $primaryTag)){
					my @tagValues = $featureObj->get_tag_values($innerTag);
					foreach my $tagToVerify (@$tagsToVerify){
						if ($tagToVerify eq $tagValues[0]){
							my %temp = exists($tagsFound{$tagToVerify}) ? $tagsFound{$tagToVerify} : ('sequence' => []);
							push($temp{'sequence'}, $featureObj->seq->seq);
							$tagsFound{$tagToVerify} = \%temp;
						}
					}
				}
			}	
		}
		$sequences{$seq->accession_number} = \%tagsFound;
	}
	return \%sequences;	
}

1;
