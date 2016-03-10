package Transformer;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Transformer);
our @EXPORT = qw(SeqIOToHashFasta SeqIOToHash organizeSequencePerORF);

sub SeqIOToHashFasta{
	my ($seqio_obj) = @_ or die "Wrong parameter number on Transformer::SeqIOToHashFasta function";
	my %sequences;
	while(my $seq = $seqio_obj->next_seq()){
		my $tempArray = ($sequences{'E1'}) ? $sequences{'E1'} : [];
		push ($tempArray, $seq->primary_seq);
		$sequences{'E1'} = $tempArray;
	}
	return \%sequences;
}

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
								unless (exists $tagsFound{$tagToVerify}) {
									$tagsFound{$tagToVerify}{'dnaSequence'} = $featureObj->seq->seq;
									my @translationArray = $featureObj->get_tag_values('translation');
									$tagsFound{$tagToVerify}{'aminoSequence'} = $translationArray[0]; 
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

sub organizeSequencesPerORF{
	my ($sequencesHash) = @_ or die "Wrong prameter on Alignment::organizeSequencesPerORF";
	my %sequences = ();
	foreach my $sequence1 (keys $sequencesHash){
		foreach my $orf (keys $sequencesHash->{$sequence1}){
			my $seq1 = $sequencesHash->{$sequence1}->{$orf}->{'dnaSequence'};
			my $tempArray = ($sequences{$orf}) ? $sequences{$orf} : [];
			push ($tempArray, Bio::PrimarySeq->new(-seq => $seq1, -id => $sequence1 . $orf));
			$sequences{$orf} = $tempArray;
		}
	}
	return \%sequences;
}

1;
