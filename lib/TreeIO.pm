package TreeIO;

use strict;
use warnings;
use Exporter;
use Data::Dumper;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(TreeIO);
our @EXPORT = qw(writeCladograms);

sub writeCladograms{
	my ($trees, $filename) = @_ or die "Wrong parameters number in createCladogram";
	my ($name, $ext) = split(/\./, $filename, 2);
	foreach my $key (keys $trees){
		my $cladogram = Bio::Tree::Draw::Cladogram->new(-bootstrap => 1, -compact => 0, -tree => $trees->{$key});
		$cladogram->print(-file => $name . "_" . $key . '.' . $ext);
	}
}

1;
