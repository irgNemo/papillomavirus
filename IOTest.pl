#!/usr/bin/perl 

use lib 'lib';
use SequenceIO;
use Bio::SeqIO;
use Transformer;
use Alignment;
use Bio::TreeIO;
use Bio::Tree::Draw::Cladogram;
use TreeIO;
use AlignmentIO;
use Data::Dumper;
use Utils;
use SequenceClustering;

my @labels = ("one", "two", "tree");
my @orfs = ('E1');#,'E2','L1','L2','E6', 'E7');
my $seqio_obj = SequenceIO::readSequence('data/sequence.gbk');
my $orfCount = Transformer::SeqIOToHash($seqio_obj, \@orfs, 'CDS', 'gene,product,note');
my $sequencesPerORF = Transformer::organizeSequencesPerORF($orfCount);

foreach my $key (keys $sequencesPerORF){
	my $seqsArrays = $sequencesPerORF->{$key};
	my $oldClusters = createRandomSets($seqsArrays, \@labels);
	my @alignmentParams = ('ktuple' => 3, 'matrix' => 'BLOSUM', 'output' => 'gcg');
	my $alignments = Alignment::clustalWAlignments($oldClusters, @alignmentParams);
	my $newClusters = SequenceClustering::computeClusters($alignments, $seqsArrays); 
	while(!areEqualSets($oldClusters, $newClusters)){
		$oldClusters = $newClusters;
		$alignments = Alignment::clustalWAlignments($oldClusters, @alignmentParams);
		$newClusters = SequenceClustering::computeClusters($alignments, $seqsArrays); 
		print "\n-- Old Cluster ---\n" . Dumper($oldClusters) . "\n-----\n";
		print "\n--- Clustered ------\n" . Dumper($newClusters) . "\n-------\n";
	}

	#print "\n --------EqualSets-------\n" . areEqualSets($randomClusters, $newClusters) . "\n----------\n";
}

sub areEqualSets{ # TODO Utilizar una algoritmo de busqueda mas eficiente
	my ($first, $second) = @_ or die "Wrong number of parameters in compareHashKeyByKey function";
	my $isFound;
	foreach my $key (keys $first){
		my $firstArray = $first->{$key};
		my $secondArray	= $second->{$key};
		print "\n------------- $key ------------\n";
		foreach my $firstSequence (@{$firstArray}){
		$isFound = 0;
			foreach my $secondSequence (@{$secondArray}){
				$isFound = $isFound || ($firstSequence->display_id eq $secondSequence->display_id) ? 1 : 0;
				print $firstSequence->display_id . " - " . $secondSequence->display_id . " - " . $isFound . "\n"; 
			}
			if (!$isFound){
				return 0;		
			} 
		}
	}
	return 1;
}

#my $alignmentsAndTrees = Alignment::clustalWAlignmentsAndTrees($sequencesPerORF, 3, 'BLOSUM');
#my @alignmentParams = ('ktuple' => 2, 'matrix' => 'BLOSUM', 'output' => 'gcg');
#my $alignments = Alignment::clustalWAlignments($sequencesPerORF, @alignmentParams);
#print Dumper($alignments);
#AlignmentIO::writeAlignments($alignments, 'msf', 'hpv.gcg');
#my @treeParams = ('outputtree' => 'nj', 'tossgaps' => 1, 'kimura' => 1, 'seed' => 121, 'bootlabels' => 'nodes', 'quiet' => 1);
#my $trees = Alignment::clustalWTrees($alignments, @treeParams);
#TreeIO::writeCladograms($trees, 'hpv.eps');
#print Dumper($alignmentsAndTrees);
#tagStatisticsReport($orfCount, \@orfs);

sub tagStatisticsReport(){
	my ($tagCount, $tagsExpected) = @_ or die "Wrong parameters number in tagStatisticReport";
	my $fila = "";
	my $header = ",";
	my $i = 0;
	foreach my $sequence (keys %$tagCount){
		$fila .= $sequence . ",";
		foreach my $tag (@$tagsExpected){
			$fila .= (exists($tagCount->{$sequence}->{$tag}->{'dnaSequence'})) ? '1' : '0';
			$fila .= ",";
			if ($i == 1){
				$header .= $tag . ",";
			}
		}
		$fila .= "\n";
		$i++;
	}
	print "$header\n $fila\n";
}

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
