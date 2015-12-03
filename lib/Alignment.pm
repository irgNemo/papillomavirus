package Alignment;

BEGIN {use Config; use Cwd; my $document_root = getcwd . "/lib/aligners/clustalw/$Config{osname}"; $ENV{CLUSTALDIR} = $document_root;}

use strict;
use warnings;
use Exporter;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::AlignIO;
use Data::Dumper;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Alignment);
our @EXPORT = qw(clustalWPairwiseAlignmentVectors clustalWAlignments clustalWTrees clustalWAlignmentsAndTrees);

sub clustalWPairwiseAlignmentVectors{
	my ($sequencesHash, @params) = @_ or die "Wrong prameter on Alignment::scoreAlignmentVectors";
	my %scoreAlignmentHash;
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
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

sub clustalWAlignments{
	my ($sequencesHash, @params) = @_ or die "Wrong prameter on Alignment::multipleAlignmentPerORF";
	
	foreach my $key (keys $sequencesHash){
		die "Less than 2 sequences in a multiple sequence alignment is not allowed" if ($#{$sequencesHash->{$key}} + 1) < 2;		
	}

	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
	my %alignments;
	foreach my $orf (keys $sequencesHash){
		$alignments{$orf} = $factory->align($sequencesHash->{$orf});
	}
	return \%alignments;
}

sub clustalWTrees{
	my ($alignments, @params) = @_ or die "Wrong prameter on Alignment::multipleAlignmentPerORF";
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
	my %trees;
	foreach my $orf (keys $alignments){
		$trees{$orf} = $factory->tree($alignments->{$orf});
	}
	return \%trees;
}

sub clustalWAlignmentsAndTrees{
	my ($sequencesHash, @params) = @_ or die "Wrong prameter on Alignment::multipleAlignmentPerORF";
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
	my %alignmentAndTree;
	foreach my $orf (keys $sequencesHash){
		my ($aln, $tree) = $factory->run($sequencesHash->{$orf});
		$alignmentAndTree{$orf}{'alignment'} = $aln;
		$alignmentAndTree{$orf}{'tree'} = $tree;
	}
	return \%alignmentAndTree;
}


1;
