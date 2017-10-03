package JackknifeBenchmarker;

use strict;
use FileHandle;

use POSIX qw(ceil floor);
use List::Util qw(sum);

my $POSITIVE = "positive";
my $NEGATIVE = "negative";

our @ISA = qw(Benchmarker);

sub new {
    my ($class, $model) = @_;
    my $self = $class->SUPER::new($model);

    bless ($self, $class);
   
    return $self;
}

#Initial policy for test set / training set values.
#Test set positives:  user-defined fraction of total positives (rounded down)
#Training set positives: all positives not in the test set
#Training set negatives: 1:1 ratio with training set positives
#Test set negatives: all remaining negatives in the test set (although we define the count explicitly beforehand and sample that many randomly for test set before creating negative training set)
sub createTestSet{

    my ($self) = @_;

    $self->writeLog("Creating test set");
    my $startTime = time();

    my $model = $self->getModel();

    #use user-defined jackknife ratio and predetermined counts of peptides of each type to determine final training, test set counts
    my $positives = $model->getPeptides($POSITIVE);
    my $negatives = $model->getPeptides($NEGATIVE);

    my $positiveCount = $model->getPeptideCount($POSITIVE);
    my $negativeCount = $model->getPeptideCount($NEGATIVE);

    my $testSetPercentage = $model->getParam("test_set_percentage");

    my $testSetPositiveCount = floor(($positiveCount * 1.0) * ($testSetPercentage * 1.0));
    my $trainingSetPositiveCount = $positiveCount - $testSetPositiveCount;
    
    my $trainingSetNegativeCount = $trainingSetPositiveCount;
    my $testSetNegativeCount = $negativeCount - $trainingSetNegativeCount;

    if ($testSetPositiveCount < 1 || $testSetNegativeCount < 1){  #this is validated on user input so should never happen
	print STDERR "JackknifeBenchmarker ERROR:  there are no peptides of a particular type in the test set (test set positive count: $testSetPositiveCount; negative count: $testSetNegativeCount\n";
	exit(1);
    }
    
    print STDERR "creating test set. Positive count (total): $positiveCount Negative count (total): $negativeCount.\nTraining set: $trainingSetPositiveCount positives and negatives.\nTest set percentage: $testSetPercentage\nPositive Count (test set): $testSetPositiveCount Negative Count (test set): $testSetNegativeCount   \n";

    #get random peptides for test set (return format $hash->{seqId}->{startPosition} = 1)
    my $testSetPositives = $self->createRandomPeptideList($positives, $testSetPositiveCount, 0);
    my $testSetNegatives = $self->createRandomPeptideList($negatives, $testSetNegativeCount, 0);
    
    $self->{TestSetPositiveCount} = $testSetPositiveCount;
    $self->{TestSetNegativeCount} = $testSetNegativeCount;

    #merge two test sets into one (format $hash->{seqId}->{startPosition} = $TYPE)
    my $finalTestSet;
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetPositives, $POSITIVE);
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetNegatives, $NEGATIVE);

    $self->{TestSet} = $finalTestSet;

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeTestSet");
}

sub createTrainingSet{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Creating Training Set");

    my $testSet = $self->{TestSet};    
    my $model = $self->getModel();
    
    my $trainingSet;  #training set is of form $trainingSet->{modbaseSeqId}->{startPosition} = $TYPE

    #get all positives not in test set
    my $positives = $model->getPeptides($POSITIVE);
    my ($positiveTrainingSet, $positiveTrainingSetCount) = $self->getRemainingPeptides($positives, $testSet);

    #get all negatives not in test set (or homologs of anything in test set)
    my $negatives = $model->getPeptides($NEGATIVE);
    my $proteinPeptideMap = $self->{PeptideProteinMap};
    my ($negativeTrainingSet, $negativeTrainingSetCount) = $self->getRemainingPeptides($negatives, $testSet, $proteinPeptideMap);

    $self->writeLog("finished creating training set.  have $positiveTrainingSetCount positives and $negativeTrainingSetCount negatives");

    #merge two training sets into one (format $hash->{seqId}->{startPosition} = $TYPE)
    $trainingSet = $self->addToPeptideSet($trainingSet, $positiveTrainingSet, $POSITIVE);
    $trainingSet = $self->addToPeptideSet($trainingSet, $negativeTrainingSet, $NEGATIVE);

    #Sanity check
    my $totalPositiveCount = $model->getPeptideCount($POSITIVE);
    my $totalNegativeCount = $model->getPeptideCount($NEGATIVE);

    my $positiveTestSetCount = $self->{TestSetPositiveCount};
    my $negativeTestSetCount = $self->{TestSetNegativeCount};

    unless (($positiveTrainingSetCount + $positiveTestSetCount == $totalPositiveCount) && ($negativeTrainingSetCount + $negativeTestSetCount == $totalNegativeCount)){
	print STDERR "ERROR: training and test set counts do not add up\n";
	print STDERR "positive training $positiveTrainingSetCount positive test $positiveTestSetCount total $totalPositiveCount\n";
	print STDERR "negative training $negativeTrainingSetCount negative test $negativeTestSetCount total $totalNegativeCount\n";  #TODO -- decide how to handle this
    }
   
    $self->{TrainingSet} = $trainingSet;
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeTrainingSet");

}


sub getRemainingPeptides{

    #This method allows us to easily implement the policy above for creating test and training sets.  
    #It is specific to that policy so if we want to make it more flexible we will need a separate workflow
    my ($self, $poolPeptides, $excludePeptides) = @_;

    my $count = 0;
    my $finalPeptides;
    foreach my $seqId (keys %$poolPeptides){
	my $startPositions = $poolPeptides->{$seqId};
	
	foreach my $startPosition (keys %$startPositions){

	    unless ($excludePeptides->{$seqId}->{$startPosition}){   #check if this peptide is in exclude
		$self->writeLog("adding sequence $seqId start position $startPosition to remaining peptides");
		$finalPeptides->{$seqId}->{$startPosition} = 1;
		$count++;
	    }
	}
    }
    return $finalPeptides;
}


sub getRemainingNonHomologousPeptides{
    
    my ($self, $poolPeptides, $excludePeptides, $proteinPeptideMap) = @_;
    my $count = 0;
    my $finalPeptides;
    
    foreach my $seqId (keys %$poolPeptides){
	my $homologs = $proteinPeptideMap->{$seqId};
	my $isHomologous;
	foreach my $homologousPeptide (keys %$homologs){
	    my ($sequenceId, $startPosition) = split('\_', $homologousPeptide);
	    $isHomologous = 1 if ($excludePeptides->{$sequenceId});
	}
	next if $isHomologous;

        my $startPositions = $poolPeptides->{$seqId};

        foreach my $startPosition (keys %$startPositions){

            unless ($excludePeptides->{$seqId}->{$startPosition}){   #check if this peptide is in exclude                                                                                                                     
		
		$self->writeLog("adding sequence $seqId start position $startPosition to remaining peptides");
		$finalPeptides->{$seqId}->{$startPosition} = 1;
		$count++;
		
	    }
        }
    }
    return ($finalPeptides, $count);


}

sub processResults {

    my ($self, $resultFh) = @_;

    #Write results of this iteration to file; will be picked up in server postprocessing

    my $startTime = time();
    $self->writeLog("Processing testing results");
    
    #use model to parse model-specific output.  SVM Model reads from file where each line is only score, so we have saved TestSetIdList to map the score to the peptide
    #TODO -- if have different model, this will probably have to be more generic
    my $testSetIdList = $self->{TestSetIdList};
    my $model = $self->getModel();
    my $resultSet = $model->parseResults($testSetIdList);

    my $fpCount = 0;
    my $tpCount = 0;
    my $foundCriticalEvalue = 0;

    my $truePositiveCount = $self->{TestSetPositiveCount};
    my $trueNegativeCount = $self->{TestSetNegativeCount};

    my $criticalEvalue;
    my $criticalTpRate;
    my $criticalFpRate;

    #each key in result id is test set identifier of form 'seq_id start_position classification'; all we need is classification to determine $fpCount and $tpCount
    foreach my $resultId (sort ({$resultSet->{$b}->{score} <=> $resultSet->{$a}->{score}} keys %$resultSet)){

	my $score = $resultSet->{$resultId}->{score};
		    
	if ($resultId =~ /$POSITIVE/){   
	    $tpCount++;
	    print $resultFh $tpCount . "\t" . $fpCount . "\t" . $score . "\n";
	}
	if ($resultId =~ /$NEGATIVE/){
	    $fpCount++;
	}

	my $tpRate = ($tpCount * 1.0) / ($truePositiveCount * 1.0);
	my $fpRate = ($fpCount * 1.0) / ($trueNegativeCount * 1.0);


	if (($tpRate >= (1.0 - $fpRate)) && ($foundCriticalEvalue == 0)){
	    
	    $foundCriticalEvalue = 1;	    
	    $criticalEvalue = $score;
	    $criticalTpRate = $tpRate;
	    $criticalFpRate = $fpRate;
	}
    }
    print $resultFh "$criticalFpRate\t$criticalTpRate\t$criticalEvalue\t$fpCount\n";  #fpCount here should be the total number of negatives in the test set

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "processResults");
}






sub createTestSet_overlap{

    my ($self) = @_;

    $self->writeLog("Creating test set");
    my $startTime = time();

    my $model = $self->getModel();

    $self->writeLog("***** Creating Test Set ******");

    my $peptideProteinMap = $self->loadPeptideProteinMap($model->getParam("peptide_protein_map_file"));
    $self->{PeptideProteinMap} = $peptideProteinMap;


    #use user-defined jackknife ratio and predetermined counts of peptides of each type to determine final training, test set counts

    my $positives = $model->getPeptides($POSITIVE);
    my $negatives = $model->getPeptides($NEGATIVE);

    my $positiveCount = $model->getPeptideCount($POSITIVE);
    my $negativeCount = $model->getPeptideCount($NEGATIVE);

    my $testSetPercentage = $model->getParam("test_set_percentage");

    my $testSetPositiveCount = floor(($positiveCount * 1.0) * ($testSetPercentage * 1.0));
    my $trainingSetPositiveCount = $positiveCount - $testSetPositiveCount;
    
    my $trainingSetNegativeCount = $trainingSetPositiveCount;
    my $testSetNegativeCount = $negativeCount - $trainingSetNegativeCount;

    if ($testSetPositiveCount < 1 || $testSetNegativeCount < 1){  #this is validated on user input so should never happen
	print STDERR "JackknifeBenchmarker ERROR:  there are no peptides of a particular type in the test set (test set positive count: $testSetPositiveCount; negative count: $testSetNegativeCount\n";
	exit(1);
    }
    
    print STDERR "creating test set. Positive count (total): $positiveCount Negative count (total): $negativeCount.\nTraining set: $trainingSetPositiveCount positives and negatives.\nTest set percentage: $testSetPercentage\nPositive Count (test set): $testSetPositiveCount Negative Count (test set): $testSetNegativeCount   \n";

    my $testSetPositives = $self->getPositivePeptideTestSet($positives, $peptideProteinMap, $testSetPositiveCount);

    


    my $negativeTestSetPool = $self->makeNegativeTestSetPool($testSetPositives, $peptideProteinMap);
    my $testSetNegatives = $self->createRandomPeptideList($negativeTestSetPool, $testSetNegativeCount, 0);

    $self->{TestSetPositiveCount} = $testSetPositiveCount;
    $self->{TestSetNegativeCount} = $testSetNegativeCount;

    #merge two test sets into one (format $hash->{seqId}->{startPosition} = $TYPE)
    my $finalTestSet;
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetPositives, $POSITIVE);
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetNegatives, $NEGATIVE);

    $self->{TestSet} = $finalTestSet;

    $self->writeLog("final test set peptide list");
    foreach my $sequenceId (keys %$finalTestSet){
	my $peptides = $finalTestSet->{$sequenceId};
	foreach my $peptide (keys %$peptides){
	    $self->writeLog("$sequenceId\t$peptide");
	    
	}
    }

    my $tempPositiveCount = $self->countPeptides($testSetPositives);
    my $tempNegativeCount = $self->countPeptides($testSetNegatives);
    print STDERR "final final test set counts: $tempPositiveCount positives and $tempNegativeCount negatives\n";
    
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeTestSet");

}




sub createTestAndTrainingSets_overlap{

    my ($self) = @_;

    $self->writeLog("Creating test set");
    my $startTime = time();

    my $model = $self->getModel();

    $self->writeLog("***** Creating Test Set ******");

    my $peptideProteinMap = $self->loadPeptideProteinMap($model->getParam("peptide_protein_map_file"));
    $self->{PeptideProteinMap} = $peptideProteinMap;

    #use user-defined jackknife ratio and predetermined counts of peptides of each type to determine final training, test set counts

    my $positives = $model->getPeptides($POSITIVE);
    my $negatives = $model->getPeptides($NEGATIVE);

    my $positiveCount = $model->getPeptideCount($POSITIVE);
    my $negativeCount = $model->getPeptideCount($NEGATIVE);

    my $testSetPercentage = $model->getParam("test_set_percentage");

    my $testSetPositiveCount = floor(($positiveCount * 1.0) * ($testSetPercentage * 1.0));
    my $trainingSetPositiveCount = $positiveCount - $testSetPositiveCount;
    
    my $trainingSetNegativeCount = $trainingSetPositiveCount;
    my $testSetNegativeCount = $negativeCount - $trainingSetNegativeCount;

    if ($testSetPositiveCount < 1 || $testSetNegativeCount < 1){  #this is validated on user input so should never happen
	print STDERR "JackknifeBenchmarker ERROR:  there are no peptides of a particular type in the test set (test set positive count: $testSetPositiveCount; negative count: $testSetNegativeCount\n";
	exit(1);
    }
    
    print STDERR "creating test set. Positive count (total): $positiveCount Negative count (total): $negativeCount.\nTraining set: $trainingSetPositiveCount positives and negatives.\nTest set percentage: $testSetPercentage\nPositive Count (test set): $testSetPositiveCount Negative Count (test set): $testSetNegativeCount   \n";

    $self->writeLog("creating test set positives");

    my $testSetPositives = $self->getHomologousSet($positives, $peptideProteinMap, $testSetPositiveCount, $POSITIVE);

    $self->writeLog("creating training set positives");

    my $trainingSetPositives = $self->getRemainingPeptides($positives, $testSetPositives);

    $self->writeLog("creating training set negatives");

    my $trainingSetNegatives = $self->getHomologousSet($negatives, $peptideProteinMap, $trainingSetNegativeCount, $NEGATIVE);

    $self->writeLog("creating test set negatives");

    my $testSetNegatives = $self->getRemainingPeptides($negatives, $trainingSetNegatives);


    $self->runConsistencyCheck($trainingSetPositives, $testSetPositives, $peptideProteinMap);
    $self->runConsistencyCheck($trainingSetNegatives, $testSetNegatives, $peptideProteinMap);

    $self->{TestSetPositiveCount} = $testSetPositiveCount;
    $self->{TestSetNegativeCount} = $testSetNegativeCount;

    #merge two test sets into one (format $hash->{seqId}->{startPosition} = $TYPE)
    my $finalTestSet;
    $self->writeLog("adding test set positives to final test set");
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetPositives, $POSITIVE);

    $self->writeLog("adding test set negatives to final test set");
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetNegatives, $NEGATIVE);

    $self->{TestSet} = $finalTestSet;




    $self->writeLog("final test set peptide list");
    foreach my $sequenceId (keys %$finalTestSet){
	my $peptides = $finalTestSet->{$sequenceId};
	foreach my $peptide (keys %$peptides){
	    $self->writeLog("$sequenceId\t$peptide");
	    
	}
    }

    my $tempPositiveCount = $self->countPeptides($testSetPositives);
    my $trainingSetNegativeCount = $self->countPeptides($testSetNegatives);
    print STDERR "final final test set counts: $tempPositiveCount positives and $trainingSetNegativeCount negatives\n";


    my $positiveTrainingSetCount = $self->countPeptides($trainingSetPositives);
    my $negativeTrainingSetCount = $self->countPeptides($trainingSetNegatives);
    $self->writeLog("final final training set counts: $positiveTrainingSetCount positives and $negativeTrainingSetCount negatives\n");



    #merge two training sets into one (format $hash->{seqId}->{startPosition} = $TYPE)
    my $trainingSet;

    $self->writeLog("adding training set positives to final training set");
    $trainingSet = $self->addToPeptideSet($trainingSet, $trainingSetPositives, $POSITIVE);

    $self->writeLog("adding training set negatives to final training set");
    $trainingSet = $self->addToPeptideSet($trainingSet, $trainingSetNegatives, $NEGATIVE);

    #Sanity check
    my $totalPositiveCount = $model->getPeptideCount($POSITIVE);
    my $totalNegativeCount = $model->getPeptideCount($NEGATIVE);

    my $positiveTestSetCount = $self->{TestSetPositiveCount};
    my $negativeTestSetCount = $self->{TestSetNegativeCount};

    unless (($positiveTrainingSetCount + $positiveTestSetCount == $totalPositiveCount) && ($negativeTrainingSetCount + $negativeTestSetCount == $totalNegativeCount)){
	print STDERR "ERROR: training and test set counts do not add up\n";
	print STDERR "positive training $positiveTrainingSetCount positive test $positiveTestSetCount total $totalPositiveCount\n";
	print STDERR "negative training $negativeTrainingSetCount negative test $negativeTestSetCount total $totalNegativeCount\n";  #TODO -- decide how to handle this
    }
   
    $self->{TrainingSet} = $trainingSet;
    
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeTestSet");

}


#test set positives
#$testSetPositives->{$sequenceId}->{$startPosition} = 1

sub makeNegativeTestSetPool{

    my ($self, $testSetPositives, $peptideProteinMap) = @_;
    $self->writeLog("making negative test set pool");
    my $negativePool;
    foreach my $sequenceId (keys %$testSetPositives){
	
	my $homologousPeptides = $peptideProteinMap->{$sequenceId};

	foreach my $homologousPeptide (keys %$homologousPeptides){
	    my ($homologSequenceId, $startPosition) = split('\_', $homologousPeptide);
	    my $type = $homologousPeptides->{$homologousPeptide};
	    if ($type eq $NEGATIVE){
		$self->writeLog("adding $homologSequenceId $startPosition to negative peptide pool");
		$negativePool->{$homologSequenceId}->{$startPosition} = 1;
	    }
	}
    }
    return $negativePool;
}

sub countPeptides{
    my ($self, $peptides)= @_;
    my $count = 0;
    foreach my $seq (keys %$peptides){
	my $peptideIds = $peptides->{$seq};
	foreach my $peptideId (keys %$peptideIds){
	    $count++;
	}
    }
    return $count;
}



sub loadPeptideProteinMap{

    my ($self, $peptideProteinMapFile) = @_; 

    my $fh = FileHandle->new("<" . $peptideProteinMapFile) || die "could not open peptide protein map file $peptideProteinMapFile\n";

    my $peptideProteinMap;

    while (<$fh>){
	chomp;
	my $line = $_;
	my @cols = split('\s', $line);
	my $sequence = $cols[0];
	for (my $i = 1; $i < scalar(@cols); $i++){
	    my $nextHomologousPeptide = $cols[$i];
	    next if ($nextHomologousPeptide =~ /mismatch/);
	    my ($homologSequenceId, $startPosition, $peptideSequence, $type) = split('\_', $nextHomologousPeptide);
	    my $homologPeptideIdentifier = $homologSequenceId . "_" . $startPosition;
	    $self->writeLog("adding next homolog $homologPeptideIdentifier for sequence $sequence");
	    $peptideProteinMap->{$sequence}->{$homologPeptideIdentifier} = $type;
	}
    }
    return $peptideProteinMap;

    #proteinPeptideMap->{$sequenceId}->{$sequenceId_$startPosition} = $type

}

sub getHomologousSet{

    my ($self, $input, $peptideProteinMap, $count, $type) = @_;

    my @peptideList = $self->createPeptideIdentifierList($input, 0);
    
    my $peptideToIndexMap = $self->createPeptideToIndexMap(@peptideList);

    my $inputSize = scalar(@peptideList);
    $self->writeLog("Loading $type homologous peptides, count $count");
    my $randomSeqs;
    my $sampledEntries;
    my $sampledEntryCount = 0;
    my $sampleCount = 0;
    while ($sampledEntryCount < $count){  
	
	#check count to avoid infinite loop
	$sampleCount++;
	if ($sampleCount > 1000000){
	    $self->writeLog("bailing out of random peptide list because sampled too much, infinite loop (count $sampleCount)\n");
	    #this avoids an infinite loop on the cluster if there are not enough peptides left to sample because of internal error.  
	    #TODO -- make nice log message that can easily be found in a search.  Considered sending email but that would probably be too complicated -- no place to write errors in output file in training mode
	    last;
	}

	#sample next entry without replacement
	my $nextSampledEntryIndex = int(rand($inputSize - 1));
	$self->writeLog("next sampled entry index: $nextSampledEntryIndex (count $sampleCount) have $sampledEntryCount in hand");
	next if ($sampledEntries->{$nextSampledEntryIndex});   #if we already have this one, start over
	my $homologousPeptides = $self->getHomologousPeptides($peptideProteinMap, $peptideList[$nextSampledEntryIndex], $type);
	
	my $homologCount;
	foreach my $peptide(@$homologousPeptides){
	    my $homologousIndex = $peptideToIndexMap->{$peptide};
	    $homologousIndex -= 1;
	    $self->writeLog("checking if homologous index $homologousIndex is present for count");
	    unless ($sampledEntries->{$homologousIndex} == 1){  #if we already put it in sampled entry count, don't include it in the homolog count here
		$self->writeLog("it was not present, incrementing count");
		$homologCount++;
	    }
	}
	

	$self->writeLog("have $homologCount homologs for this entry");
	if ($homologCount  + $sampledEntryCount + 1  > $count){  #if $@homologousPeptides == 0, and $count == $sampledEntryCount, we're on the last one and are still below threshold.
	    $self->writeLog("this puts us over the limit");
	    next;
	}  
							      	
	$sampledEntries->{$nextSampledEntryIndex} = 1;  

	my $randomEntry = $peptideList[$nextSampledEntryIndex];
	$self->writeLog("adding test set index $nextSampledEntryIndex ($randomEntry) to sampled entries");
    	

	#find peptide identifier at this index, add to list
	my ($sequenceId, $startPosition) = split('\_', $randomEntry);
	$randomSeqs->{$sequenceId}->{$startPosition} = 1;  
	$sampledEntryCount++;
	
	foreach my $homologousPeptide (@$homologousPeptides){
	    my $homologousIndex = $peptideToIndexMap->{$homologousPeptide};
	    die "did not find index for homologous peptide $homologousPeptide\n" unless $homologousIndex;
	    $homologousIndex -= 1;
	    unless ($sampledEntries->{$homologousIndex}){  #if we already put it in sampled entry count, don't include it again here
		$self->writeLog("adding test set index $homologousIndex ($homologousPeptide) to sampled entries");
		$self->writeLog("this is homologous");
		$sampledEntries->{$homologousIndex} = 1;
		my ($sequenceId, $startPosition) = split('\_', $homologousPeptide);
		$randomSeqs->{$sequenceId}->{$startPosition} = 1;  
		$sampledEntryCount++;
	    }
	}
    }

    return $randomSeqs;
}

sub getHomologousPeptides{
    my ($self, $peptideProteinMap, $peptideId, $type) = @_;
    my ($sequenceId, $startPosition) = split('\_', $peptideId);

    my $otherPeptides = $peptideProteinMap->{$sequenceId};
    my @finalPeptides;
    my $peptideCount = 0;
    foreach my $otherPeptide (keys %$otherPeptides){
	next if ($otherPeptide eq $peptideId);  #skip input peptide, add all others to list
	if ($otherPeptides->{$otherPeptide} eq $type){
	    push (@finalPeptides, $otherPeptide);
	}
    }
    return \@finalPeptides;
}


sub createPeptideToIndexMap{

    my ($self, @peptideList) = @_;

    my $map;
    my $counter = 1;
    foreach my $peptide(@peptideList){
	$map->{$peptide} = $counter;

	$counter++;

    }
    return $map;

}


sub runConsistencyCheck{

    my ($self, $trainingSet, $testSet, $peptideProteinMap) = @_;

    foreach my $sequenceId (keys %$trainingSet){
	my $mappedPeptides = $peptideProteinMap->{$sequenceId};
	foreach my $mappedPeptide (keys %$mappedPeptides){
	    my ($mappedSequenceId, $startPosition) = split('\_', $mappedPeptide);
	    if ($testSet->{$mappedSequenceId}){
		$self->writeLog("ERROR: found homolog sequence $mappedSequenceId in test set.  This is a homolog of $sequenceId which was in the training set");
	    }
	}
    }
}

sub createTrainingSet_overlap{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Creating Training Set");

    my $testSet = $self->{TestSet};    
    my $model = $self->getModel();
    
    my $trainingSet;  #training set is of form $trainingSet->{modbaseSeqId}->{startPosition} = $TYPE

    #get all positives not in test set
    my $positives = $model->getPeptides($POSITIVE);
    $self->writeLog("Creating positive training set");
    my ($positiveTrainingSet, $positiveTrainingSetCount) = $self->getRemainingPeptides($positives, $testSet);

    #get all negatives not in test set
    my $negatives = $model->getPeptides($NEGATIVE);
    $self->writeLog("Creating negative training set");
    my $peptideProteinMap = $self->{PeptideProteinMap};
    my ($negativeTrainingSet, $negativeTrainingSetCount) = $self->getRemainingNonHomologousPeptides($negatives, $testSet, $peptideProteinMap);



    $self->runConsistencyCheck($positiveTrainingSet, $testSet, $peptideProteinMap);
    $self->runConsistencyCheck($negativeTrainingSet, $testSet, $peptideProteinMap);

    
    my $tempPositiveCount = $self->countPeptides($positiveTrainingSet);
    my $tempNegativeCount = $self->countPeptides($negativeTrainingSet);
    $self->writeLog("final final training set counts: $tempPositiveCount positives and $tempNegativeCount negatives\n");



    #merge two training sets into one (format $hash->{seqId}->{startPosition} = $TYPE)
    $trainingSet = $self->addToPeptideSet($trainingSet, $positiveTrainingSet, $POSITIVE);
    $trainingSet = $self->addToPeptideSet($trainingSet, $negativeTrainingSet, $NEGATIVE);

    #Sanity check
    my $totalPositiveCount = $model->getPeptideCount($POSITIVE);
    my $totalNegativeCount = $model->getPeptideCount($NEGATIVE);

    my $positiveTestSetCount = $self->{TestSetPositiveCount};
    my $negativeTestSetCount = $self->{TestSetNegativeCount};

    unless (($positiveTrainingSetCount + $positiveTestSetCount == $totalPositiveCount) && ($negativeTrainingSetCount + $negativeTestSetCount == $totalNegativeCount)){
	print STDERR "ERROR: training and test set counts do not add up\n";
	print STDERR "positive training $positiveTrainingSetCount positive test $positiveTestSetCount total $totalPositiveCount\n";
	print STDERR "negative training $negativeTrainingSetCount negative test $negativeTestSetCount total $totalNegativeCount\n";  #TODO -- decide how to handle this
    }
   
    $self->{TrainingSet} = $trainingSet;
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeTrainingSet");

}

sub getIterationCount{

    my ($self) = @_;
    return 1;
}


sub createTrainingSet_deprecated{

    #this was in place before the new policy to specify training and test set counts
    # it is more flexible than the current policy so if we change that, this will be useful
    my ($self) = @_;

    my $testSet = $self->{TestSet};    
    my $model = $self->getModel();
    
    my $negatives = $model->getPeptides($NEGATIVE);

    my $trainingSet;  #training set is of form $trainingSet->{modbaseSeqId}->{startPosition} = $TYPE

    #Add positives to training set.  This count is automatic; just everything not in the positive test set
    my $positives = $model->getPeptides($POSITIVE);
    my $positiveTrainingSetCount = 0;
        
    foreach my $positiveSequenceId (keys %$positives){
	my $startPositions = $positives->{$positiveSequenceId};
	
	foreach my $startPosition (keys %$startPositions){

	    unless ($testSet->{$positiveSequenceId}->{$startPosition}){   #check if this cleavage sequence is in the test set
		$trainingSet->{$positiveSequenceId}->{$startPosition} = $POSITIVE;
		$positiveTrainingSetCount++;
	    }
	}
    }
    print STDERR "added $positiveTrainingSetCount positives to training set\n";
    #negatives
    my $negativeTrainingSetRatio = $model->getParam("negative_training_set_ratio");
    my $negativeTrainingSetCount = $negativeTrainingSetRatio * $positiveTrainingSetCount;
    print STDERR "adding $negativeTrainingSetCount negatives to training set\n";
    my $source;
    my $trainingSetNegatives = $self->createRandomPeptideList($negatives, $negativeTrainingSetCount, $testSet, $source);
    
    $trainingSet = $self->addToPeptideSet($trainingSet, $trainingSetNegatives, $NEGATIVE);
    
    $self->{TrainingSet} = $trainingSet;
}


return 1;
