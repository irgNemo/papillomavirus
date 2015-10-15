package Alignment;

BEGIN {use Cwd; my $document_root = getcwd . "/lib/aligners/clustalw/"; $ENV{CLUSTALDIR} = $document_root; print $document_root . "\n"; }

use strict;
use warnings;
use Exporter;
use Bio::Tools::Run::Alignment::Clustalw;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Alignment);
our @EXPORT = qw(clustalwPairwiseAlignmentVectors, clustalwMultipleAlignmentPerORF);

sub clustalwPairwiseAlignmentVectors{
	my ($sequencesHash, $ktuple, $matrix) = @_ or die "Wrong prameter on Alignment::scoreAlignmentVectors";
	my %scoreAlignmentHash;
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => $matrix, -ktuple => $ktuple);
	foreach my $sequence1 (keys $sequencesHash){
		foreach my $sequence2 (keys $sequencesHash){
			my $test1 = ($scoreAlignmentHash{$sequence2}{$sequence1}) ? 1 : 0;
			if (($sequence1 ne $sequence2) and !($test1)){
				my @scoreAlignmentVector = ();
				foreach my $orf (keys $sequencesHash->{$sequence1}){
					my $seq1 = $sequencesHash->{$sequence1}->{$orf}->{'sequence'}[0];
					my $seq2 = $sequencesHash->{$sequence2}->{$orf}->{'sequence'}[0];
					my @sequences = (
						Bio::PrimarySeq->new(-seq => $seq1, -id => $sequence1 . "-" . $orf), 
						Bio::PrimarySeq->new(-seq => $seq2, -id => $sequence2 . "-" . $orf)
						);
					push (@scoreAlignmentVector, $factory->align(\@sequences)->score());
				}
				$scoreAlignmentHash{$sequence1}{$sequence2} = \@scoreAlignmentVector;
			}
		}
	}
	
	# Elimina el par key/value que este vacio
	foreach my $sequence1 (keys(%scoreAlignmentHash)){
		if ((scalar keys $scoreAlignmentHash{$sequence1}) == 0){
			delete $scoreAlignmentHash{$sequence1};
		}	
	}

	return \%scoreAlignmentHash;
}

sub clustalwMultipleAlignmentPerORF{
	my ($sequencesHash, $ktuple, $matrix) = @_ or die "Wrong prameter on Alignment::multipleAlignmentPerORF";
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => $matrix, -ktuple => $ktuple);
	my %sequences;
	foreach my $sequence1 (keys $sequencesHash){
		foreach my $orf (keys $sequencesHash->{$sequence1}){
			my $seq1 = $sequencesHash->{$sequence1}->{$orf}->{'sequence'}[0];
			$sequences{$orf} = ($sequences{$orf}) ? push ($sequences{$orf}, Bio::PrimarySeq->new(-seq => $seq1, -id => $sequence1 . "-" . $orf) : ();
		}
	}
	
	# Elimina el par key/value que este vacio
	foreach my $sequence1 (keys(%scoreAlignmentHash)){
		if ((scalar keys $scoreAlignmentHash{$sequence1}) == 0){
			delete $scoreAlignmentHash{$sequence1};
		}	
	}

	return \%scoreAlignmentHash;
}

1;
