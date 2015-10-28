package TreeIO;

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(TreeIO);
our @EXPORT = qw(createCladogram writeCladograms);

sub createCladogram{
	my ($alignmentAndTrees) = @_ or die "Wrong parameters number in createCladogram";
	my @cladograms;
	foreach my $key (keys $alignmentAndTrees){
		my $cladogram = Bio::Tree::Draw::Cladogram->new(-tree => $alignmentAndTrees->{$key}{'tree'});
		push (@cladograms, $cladogram);	
	}
	return \@cladograms;
}

sub writeCladograms{
	my ($cladograms, $filename) = @_ or die "Wrong parameters number in writeCladograms";
	my $counter = 0;
	my ($name, $ext) = split(/\./, $filename, 2);
	foreach my $cladogram (@$cladograms){
		$cladogram->print(-file => $name . "_" . $counter . '.' . $ext);	
	}
}

1;
