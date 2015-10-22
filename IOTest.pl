#!/usr/bin/perl 

use lib 'lib';
use IOSequence;
use Bio::SeqIO;
use Transformer;
#use Alignment;
use Data::Dumper;


my @orfs = ('E1','E2','L1','L2','E6', 'E7');
$seqio_obj = IOSequence::readSequence('data/sequence.gb');
$orfCount = Transformer::SeqIOToHash($seqio_obj, \@orfs, 'CDS', 'gene,product,note');
print Dumper($orfCount);
#tagStatisticsReport($orfCount, \@orfs);

#$pairwiseAlignmentVector = Alignment::clustalwPairwiseAlignmentVectors($orfCount, 3, 'BLOSUM');
#$multipleAlignmentVector = Alignment::clustalwMultipleAlignmentPerORF($orfCount, 3, 'BLOSUM');
#print Dumper($scoreAlignmentVector);

sub tagStatisticsReport(){
	my ($tagCount, $tagsExpected) = @_ or die "Wrong parameters number in tagStatisticReport";
	my $fila = "";
	my $header = ",";
	my $i = 0;
	foreach my $sequence (keys %$tagCount){
		$fila .= $sequence . ",";
		foreach my $tag (@$tagsExpected){
			$fila .= (exists($tagCount->{$sequence}->{$tag}->{'dnaSequence'})) ? '1' : '0';
			$fila .= ",";
			if ($i == 1){
				$header .= $tag . ",";
			}
		}
		$fila .= "\n";
		$i++;
	}
	print "$header\n $fila\n";
}

sub ioTest{
	my ($seqio_obj) = @_ or die "Wrong parameters number in ioTest";
	print $seqio_obj;
	while (my $seq = $seqio_obj->next_seq){
		print $seq->description, "\n";
		foreach my $featureObj ($seq->get_SeqFeatures){
			if ($featureObj->primary_tag eq "gene"){
				foreach my $tag ($featureObj->get_all_tags){
					if($tag eq "gene"){
						print $featureObj->get_tag_values($tag), "\t", $featureObj->seq->seq;
					}
				}
				print "\n";
			}	
		}	
	}
}
