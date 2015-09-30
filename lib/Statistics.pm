package Statistics;
use strict;
use warnings;
use Exporter;
use Data::Dumper;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Statistics);
our @EXPORT = qw(statisticsPerTag);

sub statisticsPerTag{
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
							my %temp = exists($tagsFound{$tagToVerify}) ? $tagsFound{$tagToVerify} : ('sequence' => [],'count' => 0);
							push($temp{'sequence'}, $featureObj->seq->seq);
							$temp{'count'} += 1;
							$tagsFound{$tagToVerify} = \%temp;
						}
					}
				}
			}	
		}
		$sequences{$seq->description} = \%tagsFound;
	}
	return \%sequences;	
}


#print $tagToVerify, " ", $tagsValues[0], " ", ($tagToVerify eq $tagsValues[0]), "\n";
1;
