package IOSequence;
use strict;
use warnings;
use Exporter; 

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(IOSequence);
our @EXPORT = qw(readSequence writeSequence);

sub readSequence{
	my ($path) = @_ or die "Wrong parameter number in readSequence function";
	my $seqio_obj = Bio::SeqIO->new(-file => $path);
	return $seqio_obj;
}

sub writeSequence{
	print "Falta la implementacion";
}

1;
