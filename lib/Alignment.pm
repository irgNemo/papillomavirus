package Alignment;

use strict;
use warnings;
BEGIN{ 
	BEGIN {use Cwd; my $document_root = getcwd; $document_root = "$document_root" . "/lib/aligners/clustalw"; $ENV{CLUSTALDIR} = $document_root }
}
use Exporter;
use Data::Dumper;
use Bio::Tools::Run::Alignment::Clustalw;


our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Alignment);
our @EXPORT = qw(similarityVectors);

sub init{
print "Hola";

}

sub similarityVectors{
	my ($sequencesHash) = @_ or die "Wrong prameter on Alignment::similarityVectors";
	foreach my $sequence1 (keys $sequencesHash){
		foreach my $sequence2 (keys $sequencesHash){
			if ($sequence1 ne $sequence2){
				foreach my $orf (keys $sequencesHash->{$sequence1}){
					if ($sequence1 ne $sequence2){
						my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => 'BLOSUM');
						my $ktuple = 3;
						my $seq1 = $sequencesHash->{$sequence1}->{$orf}->{'sequence'}[0];
						my $seq2 = $sequencesHash->{$sequence2}->{$orf}->{'sequence'}[0];
						my @sequences = (
							Bio::PrimarySeq->new(-seq => $seq1, -id => "$sequence1" . "-" . "$orf"), 
							Bio::PrimarySeq->new(-seq => $seq2, -id => "$sequence2" . "-" . "$orf")
							);
						my $aln = $factory->align(\@sequences);
					}
				}	
			}
		}
	}
}
