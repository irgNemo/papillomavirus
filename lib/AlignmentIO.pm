package AlignmentIO;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(AlignmentIO);
our @EXPORT = qw(writeAlignments readAlignments);

sub writeAlignments{
	my ($alignments, $format, $filename) = @_ or die "Wrong parameters number in writeAlignment function";
	my ($name, $ext) = split (/\./, $filename, 2);
	foreach my $key (keys $alignments){
		my $out = Bio::AlignIO->new(-file => ">$name" . "_" . $key . "." . $ext, '-format' => $format);
		$out->write_aln($alignments->{$key});	
	}	
}

sub readAlignments{
	my ($alignmentFile) = @_ or die "Wrong parameters number in readAlignment function";
	my %alignments;
	my $alignObject = Bio::AlignIO->new(-file => $alignmentFile);
	foreach my $aln ($alignObject->next_aln){
		$alignments{$aln->id} = $aln;
	}
	return \%alignments;
}

1;
