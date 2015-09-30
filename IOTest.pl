#!/usr/bin/perl 

use lib 'lib';
use IOSequence;
use Bio::SeqIO;
use Statistics;
use Data::Dumper;

my @orfs = ('E1','E2','L1','L2','N1');
$seqio_obj = IOSequence::readSequence('data/sequence.gb');
$p = Statistics::statisticsPerTag($seqio_obj, \@orfs, 'CDS', 'gene');
print Dumper($p);
#ioTest($seqio_obj);

sub ioTest{
	my ($seqio_obj) = @_ or die "Wrong parameters number in ioTest";
	print $seqio_obj;
	while (my $seq = $seqio_obj->next_seq){
		print $seq->description, "\n";
		foreach my $featureObj ($seq->get_SeqFeatures){
			if ($featureObj->primary_tag eq "gene"){
				foreach my $tag ($featureObj->get_all_tags){
					if($tag eq "gene"){
						print $featureObj->get_tag_values($tag), "\t", $featureObj->seq->seq;
					}
				}
				print "\n";
			}	
		}	
	}
}
