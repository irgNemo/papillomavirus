package Alignment;

BEGIN {use Config; use Cwd; my $document_root = getcwd . "/lib/aligners/clustalw/$Config{osname}"; $ENV{CLUSTALDIR} = $document_root;}

use strict;
use warnings;
use Exporter;
use Bio::Tools::Run::Alignment::Clustalw;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Alignment);
our @EXPORT = qw(clustalwPairwiseAlignmentVectors clustalwMultipleAlignment organizeSequencesPerORF clustalwAlignmentAndTree);

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

sub organizeSequencesPerORF{
	my ($sequencesHash) = @_ or die "Wrong prameter on Alignment::organizeSequencesPerORF";
	my %sequences = ();
	foreach my $sequence1 (keys $sequencesHash){
		foreach my $orf (keys $sequencesHash->{$sequence1}){
			my $seq1 = $sequencesHash->{$sequence1}->{$orf}->{'sequence'}[0];
			my $tempArray = ($sequences{$orf}) ? $sequences{$orf} : [];
			push ($tempArray, Bio::PrimarySeq->new(-seq => $seq1, -id => $sequence1 . "-" . $orf));
			$sequences{$orf} = $tempArray;
		}
	}
	return \%sequences;
}

sub clustalwMultipleAlignment{
	my ($sequencesHash, $ktuple, $matrix) = @_ or die "Wrong prameter on Alignment::multipleAlignmentPerORF";
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => $matrix, -ktuple => $ktuple);
	foreach my $orf (keys $sequencesHash){
		$factory->align($sequencesHash->{$orf});
	}
}

sub clustalwAlignmentAndTree{
	my ($sequencesHash, $ktuple, $matrix) = @_ or die "Wrong prameter on Alignment::multipleAlignmentPerORF";
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => $matrix, -ktuple => $ktuple);
	my %alignmentAndTree;
	foreach my $orf (keys $sequencesHash){
		my ($aln, $tree) = $factory->run($sequencesHash->{$orf});
		$alignmentAndTree{$orf}{'alignment'} = $aln;
		$alignmentAndTree{$orf}{'tree'} = $tree;
	}
	return \%alignmentAndTree;
}

1;
