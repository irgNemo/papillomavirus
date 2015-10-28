package Alignment;

BEGIN {use Config; use Cwd; my $document_root = getcwd . "/lib/aligners/clustalw/$Config{osname}"; $ENV{CLUSTALDIR} = $document_root;}

use strict;
use warnings;
use Exporter;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::AlignIO;

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
					my $seq1 = $sequencesHash->{$sequence1}->{$orf}->{'dnaSequence'};
					my $seq2 = $sequencesHash->{$sequence2}->{$orf}->{'dnaSequence'};
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

sub clustalwMultipleAlignment{
	my ($sequencesHash, $ktuple, $matrix) = @_ or die "Wrong prameter on Alignment::multipleAlignmentPerORF";
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => $matrix, -ktuple => $ktuple);
	my $alignment;
	foreach my $orf (keys $sequencesHash){
		$alignment = $factory->align($sequencesHash->{$orf});
	}
	return $alignment;
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

sub writeAlignment{
	my ($alignment,$filename) = @_ or die "Wrong parameters number in writeAlignment function";
	my $out = Bio::AlignIO->new(-file => ">$filename");
	while (my $aln = $alignment->next_aln){
		$out->write_aln($aln);	
	}
}

1;
