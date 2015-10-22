package Transformer;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Transformer);
our @EXPORT = qw(SeqIOToHash);

sub SeqIOToHash{
	my ($seqio_obj, $tagsToVerify, $primaryTag, $secondaryTags) = @_ or die "Wrong parameter number on Statistic::statisticPerTag function";
	my %sequences;
	while(my $seq = $seqio_obj->next_seq){
		my %tagsFound;
		foreach my $featureObj ($seq->get_SeqFeatures){
			foreach my $secondaryTag (split (",",$secondaryTags)){
				foreach my $innerTag ($featureObj->get_all_tags){
					if (($secondaryTag eq $innerTag) and ($featureObj->primary_tag eq $primaryTag)){
						my @tagValues = $featureObj->get_tag_values($innerTag);
						foreach my $tagToVerify (@$tagsToVerify){
							my ($tagFounded) = (($tagValues[0] =~ m/\s+($tagToVerify)\s+/) || ($tagValues[0] =~ m/^($tagToVerify)$/) || ($tagValues[0] =~ m/\s+($tagToVerify)$/) || ($tagValues[0] =~ m/^($tagToVerify)\s+/));
							if (defined($tagFounded)){
								unless (exists $tagsFound{$tagToVerify}) {															$tagsFound{$tagToVerify}{'dnaSequence'} = $featureObj->seq->seq;
								}
							}
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
