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

sub initialize{  #currently nothing
    my ($self) = @_;
}


#######################################################################################################################################
# createTestSet()
# Takes peptides loaded by model and partitions a random selection of them into a test set. 
# The number of each type (positive and negative) to put in the test set is calcluated as follows:
#
# 1. Test Set Positive count: Apply user defined test_set_percentage (i.e., jackknife fraction) to total number of positives; round down
#    Randomly select positives to put into test set
# 2. Training set positive count: All positives not in the test set go here
# 3. Training set negative count: Defined as an equal number to the number of training set positives; randomly selected into training set
# 4. Test set positive count: Take all remaining negatives and put them into the test set.
#
# While all counts are defined in this method, only the test set peptides are initialized here. They are stored in $self->{TestSet} in a hash
# of the form:
# 
# $testSet->{$modbaseSeqId}->{$startPosition} = $classification (positive/negative)
#
# RETURN NULL
#######################################################################################################################################
sub createTestSet{

    my ($self) = @_;

    $self->writeLog("Creating test set");
    my $startTime = time();

    return 0 if ($self->pipelineHasErrors());

    my $model = $self->getModel();

    #Calculate training and test set counts
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
	$model->writeError("invalid_benchmark_ratio", "createTestSet", "Test set peptide count less than 0 ($testSetPositiveCount positives, $testSetNegativeCount negatives)");
	return 0;
    }
    
    $self->writeLog("Benchmark statistics: $positiveCount total positives ($trainingSetPositiveCount training peptides and $testSetPositiveCount test peptides)");
    $self->writeLog("$negativeCount total negatives ($trainingSetNegativeCount training peptides and $testSetNegativeCount test peptides)");

    #get random peptides for test set (return format $hash->{seqId}->{startPosition} = 1)
    my $testSetPositives = $self->createRandomPeptideList($positives, $testSetPositiveCount, 0);
    my $testSetNegatives = $self->createRandomPeptideList($negatives, $testSetNegativeCount, 0);
    
    return 0 unless ($testSetPositives && $testSetNegatives);

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

#######################################################################################################################################
# createTrainingSet()
# Takes peptides loaded by model and partitions a random set of them into a training set.
# Since the test set peptides have already been loaded, this method gets everything not already in the test set and puts it in the 
# training set. Peptides are stored in $self->{TrainingSet} which is of the form
#
# $trainingSet->{$modbaseSeqId}->{$startPosition} = $classification (positive/negative)
#
# RETURN NULL
#######################################################################################################################################
sub createTrainingSet{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Creating Training Set");

    return 0 if ($self->pipelineHasErrors());

    my $testSet = $self->{TestSet};    
    my $model = $self->getModel();
    
    my $trainingSet;  

    #get all positives not in test set
    my $positives = $model->getPeptides($POSITIVE);
    my ($positiveTrainingSet, $positiveTrainingSetCount) = $self->getRemainingPeptides($positives, $testSet);
    
    #get all negatives not in test set 
    my $negatives = $model->getPeptides($NEGATIVE);
    my ($negativeTrainingSet, $negativeTrainingSetCount) = $self->getRemainingPeptides($negatives, $testSet);

    #merge two training sets into one
    $trainingSet = $self->addToPeptideSet($trainingSet, $positiveTrainingSet, $POSITIVE);
    $trainingSet = $self->addToPeptideSet($trainingSet, $negativeTrainingSet, $NEGATIVE);

    #Sanity check; make sure training and test set observed counts add up to total peptide counts
    my $totalPositiveCount = $model->getPeptideCount($POSITIVE);
    my $totalNegativeCount = $model->getPeptideCount($NEGATIVE);
    my $positiveTestSetCount = $self->{TestSetPositiveCount};
    my $negativeTestSetCount = $self->{TestSetNegativeCount};
    
    unless (($positiveTrainingSetCount + $positiveTestSetCount == $totalPositiveCount) && ($negativeTrainingSetCount + $negativeTestSetCount == $totalNegativeCount)){
	my $errorMsg = "Calculations do not add up after compiling all benchmark sets. Observed counts: $totalPositiveCount positives ($positiveTrainingSetCount training and ";
	$errorMsg .= "$positiveTestSetCount test); $totalNegativeCount negatives ($negativeTrainingSetCount training and $negativeTestSetCount test)";
	$model->writeError("invalid_benchmark_ratio", "createTrainingSet", $errorMsg);
    }

    $self->writeLog("Done creating training set");

    $self->{TrainingSet} = $trainingSet;
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeTrainingSet");

}

#######################################################################################################################################
# processResults()
# Create output for this benchmark run. Uses the model to read its output score file, then gets a list of all peptides in the test
# set, sorted by score. With this list, output one line in the result file for each positive peptide. The line is tab delimited
# and is of the form (tpCount\tfpCount\tscore). This represents the number of true positives and false positives that are at or better
# than the score. Note that the $tpCount column will increase incrementally. 
#
# The bottom line of the file includes the critical fp and tp rates, and the total number of negatives in the test set. 
#
# PARAM  $resultFh: The FileHandle to which to write results.
# RETURN $errors if any occured; else return 0
#######################################################################################################################################
sub processResults {

    my ($self, $resultFh) = @_;

    #Check if errors; if any exist, will be handled by Pipeline
    my $errors = $self->pipelineHasErrors();
    if ($errors){
	return $errors;
    }

    my $startTime = time();
    $self->writeLog("Processing Results");
    
    #Assign scores generated by the model to all peptides, using $testSetIdList (set by Benchmarker in testModel())
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

    print $resultFh "TP_count\tFP_count\tScore\n";

    #Go down through list of peptides, sorted by scores, to determine TPR, FPR at each level 
    foreach my $resultId (sort ({$resultSet->{$b}->{score} <=> $resultSet->{$a}->{score}} keys %$resultSet)){

	my $score = $resultSet->{$resultId}->{score};

	#$resultId: form of 'seq_id start_position classification'; all we need is classification to determine $fpCount and $tpCount		    
	if ($resultId =~ /$POSITIVE/){   
	    $tpCount++;
	    print $resultFh $tpCount . "\t" . $fpCount . "\t" . $score . "\n";
	}
	if ($resultId =~ /$NEGATIVE/){
	    $fpCount++;
	}

	#update current rates
	my $tpRate = ($tpCount * 1.0) / ($truePositiveCount * 1.0);
	my $fpRate = ($fpCount * 1.0) / ($trueNegativeCount * 1.0);

	#Check if we've gone past the 'critical point'
	if (($tpRate >= (1.0 - $fpRate)) && ($foundCriticalEvalue == 0)){
	    
	    $foundCriticalEvalue = 1;	    
	    $criticalEvalue = $score;
	    $criticalTpRate = $tpRate;
	    $criticalFpRate = $fpRate;
	}
    }
    print $resultFh "$criticalFpRate\t$criticalTpRate\t$criticalEvalue\t$fpCount\n";  
    #fpCount here should be the total number of negatives in the test set. Sort of hacky but couldn't figure out another way to communicate
    #that to the backend
    
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "processResults");
    return 0;
}


#######################################################################################################################################
# getRemainingPeptides
# Given two sets of peptides, return all peptides that are in one but not the other.
# 
# PARAM  $poolPeptides: The larger pool of peptides (format: $poolPeptides->{$modbaseSeqId}->{$startPosition} = $seqInfo; $seqInfo unused here)
# PARAM  $excludePeptides: The smaller pool of peptides (format $excludePeptides->{$modbaseSeqId}->{$startPosition} = $classification; $classification unused here)
# RETURN All peptides that are in $poolPeptides but not $excludePeptides. Returns a hash of the form
#        $finalPeptides->{$modbaseSequenceId}->{$startPosition} = 1
#######################################################################################################################################
sub getRemainingPeptides{

    my ($self, $poolPeptides, $excludePeptides) = @_;

    my $finalPeptides;
    my $count = 0;
    my $debugCount = 0;
    my $excludeCount = 0;
    foreach my $seqId (keys %$poolPeptides){
	my $startPositions = $poolPeptides->{$seqId};
	
	foreach my $startPosition (keys %$startPositions){
	    $debugCount++;
	    unless ($excludePeptides->{$seqId}->{$startPosition}){   #check if this peptide is in exclude
		$finalPeptides->{$seqId}->{$startPosition} = 1;
		$count++;
	    }
	    else {
		$excludeCount++;
	    }
	}
    }
    return ($finalPeptides, $count);
}

############################################################################################################################################################     
# pipelineHasErrors
# 
# Checks to make sure that no errors occured in that previous pipeline or in this one. 
#
# RETURN: The name of the error if a fatal error occured, OR 0 if everything worked normally.
############################################################################################################################################################     
sub pipelineHasErrors{
    my ($self) = @_;
    my $model = $self->getModel();

    my $errors = $model->getErrors();
    if ($errors){
	return $errors;
    }

    my $pipelineErrors = $self->getErrors();
    if ($pipelineErrors){
	return $pipelineErrors;
    }

    return 0;
}




sub setErrors{
    my ($self, $errorKeyword) = @_;
    $self->{errors} = $errorKeyword;
}

sub getErrors{
    my ($self) = @_;
    my $errors = $self->{errors};
    return $errors;
}



sub getIterationCount{

    my ($self) = @_;
    return 1;
}





###################################################################################################################################################
###################################################################################################################################################
#Everything below this was only used for the paper, and not in the web service






###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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
		$finalPeptides->{$seqId}->{$startPosition} = 1;
		$count++;
		
	    }
        }
    }
    return ($finalPeptides, $count);


}



###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################

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




###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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

###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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

###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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


###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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

###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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
###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################

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

###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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

###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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

###################################################################################################################################################
# Method created to satisfy concern of reviewers that there are no problems with training / test set peptides from homologous sequences
###################################################################################################################################################
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
