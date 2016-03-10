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

my @labels = ("one", "two", "tree", "four", "five");
my @orfs = ('E1');#,'E2','L1','L2','E6', 'E7');
my $seqio_obj = SequenceIO::readSequence('data/sequence.gbk'); # Este lee archivos gb
#my $seqio_obj = SequenceIO::readSequence('data/sequence.gbk');
my $orfCount = Transformer::SeqIOToHash($seqio_obj, \@orfs, 'CDS', 'gene,product,note');
my $sequencesPerORF = Transformer::organizeSequencesPerORF($orfCount);
my $stopCondition = 10;
my @alignmentParams = ('ktuple' => 3, 'matrix' => 'BLOSUM', 'output' => 'gcg', 'quiet' => '1');
my $clustersAndAlignments = undef;
my $maxClusteringIterations = 10; # Numero de iteraciones clustering
my $currentClusteringIteration = 0; # Actual iteracion
my $finalClusterAndAlignment = undef;
my $maxAverageScore = 0;


my $start = DateTime->now();
while($currentClusteringIteration < $maxClusteringIterations){
	print "\n----------- Current Ieration $currentClusteringIteration ------------\n";
	print "\n----------------- SequenceClustering --------------------\n";
	$clustersAndAlignments = sequenceClustering($sequencesPerORF, \@alignmentParams, $stopCondition);
	my $averageScoreAlignment = averageScoreAlignmentPerCluster($clustersAndAlignments->{"alignments"});
	if ($maxAverageScore < $averageScoreAlignment) {
		$maxAverageScore = $averageScoreAlignment;
		$finalClusterAndAlignment = $clustersAndAlignments; 
	}
	print "\n----------------- Average Score $averageScoreAlignment --------------------\n";
	$currentClusteringIteration++;
}

my $end = DateTime->now();
my $elapse = $end - $start;
print "\n------------------Cluster Final -------------------\n";
imprimirClusters($finalClusterAndAlignment->{"clusters"});
print "\n--------- Average Score:" .  averageScoreAlignmentPerCluster($finalClusterAndAlignment->{"alignments"}) . "-------\n";
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
		$oldClusters = SequenceClustering::createRandomSets($seqsArrays, \@labels);
		$alignments = Alignment::clustalWAlignments($oldClusters, @alignmentParams);
		print "\n-------- Inicio Iteracion $iteration ------------\n";
		$newClusters = SequenceClustering::computeClusters($alignments, $seqsArrays, @alignmentParams);
		$newClusters = validateSizeClusters($newClusters, 2); # Verificar que los nuevos Clusters tenga al menos dos secuencias por cada uno, si no, asignarlo de manera aleatoria a alguno
		print "\n --------- Old Clusters ------------\n";
		imprimirClusters($oldClusters);
		print "\n ---------- New Clusters --------\n";
		imprimirClusters($newClusters);
		print "\n-------- Fin Iteracion $iteration ------------\n";
		$iteration++;
		while(!areEqualSets($oldClusters, $newClusters) && $iteration < $stopCondition){
			$oldClusters = $newClusters;
			$alignments = Alignment::clustalWAlignments($oldClusters, @alignmentParams);
			print "\n-------- Inicio Iteracion $iteration ------------\n";
			$newClusters = SequenceClustering::computeClusters($alignments, $seqsArrays, @alignmentParams); 
			$newClusters = validateSizeClusters($newClusters, 2); # Verificar que los nuevos Clusters tenga al menos dos secuencias por cada uno, si no, asignarlo de manera aleatoria a alguno
			print "\n --------- Old Clusters ------------\n";
			imprimirClusters($oldClusters);
			print "\n ---------- New Clusters --------\n";
			imprimirClusters($newClusters);
			print "\n-------- Fin Iteracion $iteration ------------\n";
			$iteration++;
		}
	}
	return {"alignments" => $alignments, "clusters" => $newClusters};
}

sub validateSizeClusters{
	my ($clusters, $minNumOfElements) = @_ or die "Wrong number of parameters in validateSizeClusters function";
	my @labels = keys $clusters;
	foreach $label (@labels){
		my $size = @{$clusters->{$label}};
		if ($size < $minNumOfElements){
			my $randomLabel = $labels[int(rand(@labels))];
			while ($label eq $randomLabel) { $randomLabel = $labels[int(rand(@labels))];}
			print "\nnew cluster From $label To  $randomLabel\n";
			foreach my $element (@{$clusters->{$label}}){
				my @temp = @{$clusters->{$randomLabel}};
				push (@temp, $element);
				$clusters->{$randomLabel} = \@temp;
			}
			delete $clusters->{$label};
		}
	}
	return $clusters;
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
