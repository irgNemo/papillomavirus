#!/usr/bin/perl 

use lib 'lib';
use Utils;
use IOSequence;
use Data::Dumper;
use Transformer;
use Bio::SeqIO;



my @orfs = ('E6', 'E7', 'E1', 'E2', 'L2', 'L1');
my @list = ('NC_001526','AY262282','J04353','M12732','M74117','M62849','KC470260','HQ537751','D90400','D21208','DQ344807','U31791','U31788','AF436130','FR872717','AF092932','M73236');
my $path = "data/";

my $o=1;
my $l=0;

foreach $id (@list){
#get the gb files for each HPV in list
	#Utils::getNCBICompleteGenomebyID($id,$path);
}

foreach $id (@list){
	$l++;
	$o=1;
	$file = $path.$id.".gb";
	$seqio_obj = IOSequence::readSequence($file);
	$orfCount = Transformer::SeqIOToHash($seqio_obj, \@orfs, 'CDS', 'gene');
	foreach $g (@orfs){
		my $file = "data/gerardo_fasta/".$l."_".$o++.".fasta";
		unless(open FILE, '>'.$file) {
    		die "\nUnable to create $file\n";
		}
		print FILE ">HPV|$id|ORF|$g\n".$orfCount->{$id}{$g}{'sequence'}[0]."\n";
		close FILE;
	}
	#print Dumper($orfCount);
}


