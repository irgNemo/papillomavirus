package SequenceClustering;

use strict;
use warnings;
use Exporter;
use Data::Dumper;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(SequenceClustering);
our @EXPORT = qw(createRandomSets computeCluster);

# TODO Esta funcion debera tene el proceso de clustering de secuencias completo

#/////////////////////////////////////////////////////////////////////////////////////////
=head2 createRandomSets
 Parameters	: setToLabel An array reference of the elements to be rearrange
		  labels An string array with the labels for each group; at least a 2-dimensional array
		  isDisjoint True if the sets must be disjount, otherwise false.	
 Returns	: An array with as many reference sets as labels in the setLabels array. Disjoint or not depending on isDisjoint flag.
 Description	: Create as many sets as labels in setsLabels array.

=cut

sub createRandomSets{
	my ($elements, $labels) = @_ or die "Wrong parameters number in createRandomSets";
	my %sets;
	my $minSequencesPerCluster = 2;

	# Completamos los clusters, aleatoreamente, con el resto de los elementos.
	while (1){	
		foreach my $element (@{$elements}){
			my $randomIndex = int(rand(@{$labels}));
			my $label = $labels->[$randomIndex];
			my $temp = (exists $sets{$label}) ? $sets{$label} : [];
			push ($temp, $element);
			$sets{$label} = $temp;
		}
		
		my $biggerThanTwo = 1;
		foreach my $key (keys %sets){
			my $size =  @{$sets{$key}};
			$biggerThanTwo &&= ($size >= $minSequencesPerCluster) ? 1 : 0;
		}
		if ($biggerThanTwo) { last; }
		%sets = ();
	}
	return \%sets;
}

sub computeClusters{
	my ($alignments, $sequences, @alignmentParams) = @_ or die "Wrong parameters number in reassignLabels";
	my %clusters;
	foreach my $sequence (@{$sequences}){
		my $maxAlignment = 0; # Si vamos a utilizar score entonces hay que poner el valor menor (mas negativo) posible
		my $newLabel;
		my $sequenceToBeStored;
		foreach my $key (keys $alignments){
			# TODO Eliminar los objetos creados explicitamente
			my $consensus = Bio::PrimarySeq->new(-seq => $alignments->{$key}->consensus_string(60), -id => $alignments->{$key}->id()); # Ver si puedo obtener un objeto secuencia del alineamiento. Pasar como parametro el valor de la cota minima para determinar la secuencia consenso.
			my $alignment = (Alignment::clustalWAlignments({'pairwise' => [$sequence, $consensus]}, @alignmentParams)->{'pairwise'});
			my $alignmentPercentage = $alignment->percentage_identity; #average_percentage_identity(); # Verificar cual es el valor de identidad que deseamos utilizar
			if ($maxAlignment <= $alignmentPercentage or !defined $maxAlignment){
				$newLabel = $key;
				$maxAlignment = $alignmentPercentage;
				$sequenceToBeStored = $sequence;
			}
		}
		my $temp = (exists $clusters{$newLabel}) ? $clusters{$newLabel} : [];
		push($temp, $sequenceToBeStored);
		$clusters{$newLabel} = $temp;	
	}
	return \%clusters;
}

1;
