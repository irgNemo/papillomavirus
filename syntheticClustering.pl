#!/usr/bin/perl

use lib 'lib';
use SequenceIO;
use Bio::SeqIO;
use Data::Dumper;
use Transformer;

my $seqio_obj = SequenceIO::readSequence("data/secuencia1-percentage1.fasta");
my $seqHash = SeqIOToHashFasta($seqio_obj);
print Dumper($seqHash);
