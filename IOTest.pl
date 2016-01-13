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
use DateTime;

my @labels = ("one", "two");
my @orfs = ('E1');#,'E2','L1','L2','E6', 'E7');
my $seqio_obj = SequenceIO::readSequence('data/sequence.gbk');
my $orfCount = Transformer::SeqIOToHash($seqio_obj, \@orfs, 'CDS', 'gene,product,note');
my $sequencesPerORF = Transformer::organizeSequencesPerORF($orfCount);
my $stopCondition = 10;
my $maxClustersRuns = 1;
my $clustersRuns = 0;
my @alignmentParams = ('ktuple' => 3, 'matrix' => 'BLOSUM', 'output' => 'gcg', 'quiet' => '1');
my $start = DateTime->now();
my $clustersAndAlignments = undef;
my $maxCluster = 0;
my $finalClusterAndAlignment = undef;
while($clustersRuns < $maxClustersRuns){
	$clustersAndAlignments = sequenceClustering($sequencesPerORF, \@alignmentParams, $stopCondition);
	my $averageScoreAlignment = averageScoreAlignmentPerCluster($clustersAndAlignments->{"alignments"});
	if ($maxClusters < $averageScoreAlignment) {
		$maxClusters = $averageScoreAlignment;
		$finalClusterAndAlignment = $clustersAndAlignments; 
	}
	$clustersRuns++;
}

my $end = DateTime->now();
my $elapse = $end - $start;
print "------------------Cluster Final -------------------";
imprimitClusters($finalClustersAndAlignment->{"clusters"});
print "\n------------ Elapsed time : " . $elapse->in_units('minutes') . " min --------------\n";

sub averageScoreAlignmentPerCluster{
	my ($clustersAlignment) = @_ or die "Wrong number of parameters in averageScoreAlignmentPerCluser";
	my $sum = 0;
	my $numberOfClusters = 0;
	foreach my $key (keys $clustersAlignment){
		$sum = $sum + $clustersAlignment->{$key}->percentage_identity;
		$numberOfClusters++;
	}
	
	return 	($sum / ($numberOfClusters));
}

sub sequenceClustering{
	my ($sequencesPerORF, $alignmentParams, $stopCondition) = @_ or die "Wrong number of parameters in computeClusters";
	my $oldClusters = undef;
	my $newClusters = undef;
	my $alignments = undef;
	foreach my $key (keys $sequencesPerORF){
		print "\n-------------------- Llave $key ---------------\n";
		my $iteration = 0;
		my $seqsArrays = $sequencesPerORF->{$key};
		$oldClusters = createRandomSets($seqsArrays, \@labels);
		$alignments = Alignment::clustalWAlignments($oldClusters, @alignmentParams);
		print "\n-------- Iteracion $iteration ------------\n";
		$newClusters = SequenceClustering::computeClusters($alignments, $seqsArrays, @alignmentParams);
		print "\n --------- Old Clusters ------------\n";
		imprimirClusters($oldClusters);
		print "\n ---------- New Clusters --------\n";
		imprimirClusters($newClusters);
		print "\n-------- Iteracion $iteration ------------\n";
		$iteration++;
		while(!areEqualSets($oldClusters, $newClusters) && $iteration < $stopCondition){
			$oldClusters = $newClusters;
			$alignments = Alignment::clustalWAlignments($oldClusters, @alignmentParams);
			print "-------- Iteracion $iteration ------------\n";
			$newClusters = SequenceClustering::computeClusters($alignments, $seqsArrays, @alignmentParams); 
			print "\n --------- Old Clusters ------------\n";
			imprimirClusters($oldClusters);
			print "\n ---------- New Clusters --------\n";
			imprimirClusters($newClusters);
			print "-------- Iteracion $iteration ------------\n";
			$iteration++;
		}
	}
	return {"alignments" => $alignments, "clusters" => $newClusters};
}

sub imprimirClusters{
	my ($clusters) = @_ or die "Wrong number of parameters in imprimirClusters funtion";
	for my $key (keys $clusters){
		my $sequences = $clusters->{$key};
		print "----- Cluster $key -----\n";
		for my $sequence (@{$sequences}){
			print $sequence->display_id . "\n";
		}
	}
}

sub areEqualSets{ # TODO Utilizar una algoritmo de busqueda mas eficiente
	my ($first, $second) = @_ or die "Wrong number of parameters in compareHashKeyByKey function";
	my $isFound = undef;
	foreach my $key (keys $first){
		my $firstArray = $first->{$key};
		my $secondArray	= $second->{$key};

		if ($#{$firstArray} != $#{$secondArray}){
			return 0; 
		}
		foreach my $firstSequence (@{$firstArray}){
		$isFound = 0;
			foreach my $secondSequence (@{$secondArray}){
				$isFound = $isFound || ($firstSequence->display_id eq $secondSequence->display_id) ? 1 : 0;
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
