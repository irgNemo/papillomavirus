package Alignment;

use strict;
use warnings;
BEGIN{ 
	BEGIN {use Cwd; my $document_root = getcwd . "/lib/aligners/clustalw"; $ENV{CLUSTALDIR} = $document_root }
}
use Exporter;
use Data::Dumper;
use Bio::Tools::Run::Alignment::Clustalw;


our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Alignment);
our @EXPORT = qw(scoreAlignmentVectors);

sub scoreAlignmentVectors{
	my ($sequencesHash, $ktuple, $matrix) = @_ or die "Wrong prameter on Alignment::similarityVectors";
	my %scoreAlignmentHash;
	foreach my $sequence1 (keys $sequencesHash){
		foreach my $sequence2 (keys $sequencesHash){
			if (($sequence1 ne $sequence2) or (!(exists($scoreAlignmentHash{$sequence1}{$sequence2})) and !(exists($scoreAlignmentHash{$sequence2}{$sequence1})))){
				my @scoreAlignmentVector = ();
				foreach my $orf (keys $sequencesHash->{$sequence1}){
					my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => $matrix, -ktuple => $ktuple);
					$factory->cleanup();
					my $seq1 = $sequencesHash->{$sequence1}->{$orf}->{'sequence'}[0];
					my $seq2 = $sequencesHash->{$sequence2}->{$orf}->{'sequence'}[0];
					my @sequences = (
						Bio::PrimarySeq->new(-seq => $seq1, -id => "$sequence1" . "-" . "$orf"), 
						Bio::PrimarySeq->new(-seq => $seq2, -id => "$sequence2" . "-" . "$orf")
						);
					push (@scoreAlignmentVector, $factory->align(\@sequences)->score());
				}
				$scoreAlignmentHash{$sequence1}{$sequence2} = \@scoreAlignmentVector;
			}
		}
	}
	return \%scoreAlignmentHash;
}
