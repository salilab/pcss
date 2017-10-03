package LeaveOneOutBenchmarker;

use strict;
use FileHandle;

my $POSITIVE = "positive";
my $NEGATIVE = "negative";
my $TEST = "Test";
my $PROTEOME = "Proteome";



use POSIX qw(ceil floor);
use List::Util qw(sum);

our @ISA = qw(Benchmarker);

sub new {
    my ($class, $model) = @_;
    my $self = $class->SUPER::new($model);

    bless ($self, $class);
    $self->{PeptidePoolLoaded} = 0;
   
    return $self;
}


#######################################################################################################################################
# initialize
# All this does is call initializePeptidePool(). This method should be called AFTER the peptides have been loaded into the Model.
# However if it's not, createTestSet() will do so.
#######################################################################################################################################
sub initialize{
    my ($self) = @_;
   
    $self->initializePeptidePool(); 
}



#######################################################################################################################################
# initializePeptidePool
# Reads in training input and adds to internal lists of positives and negatives.
# Lists are stored in an internal hash of the form:
# 
# $self->{PeptidePool}->{$POSITIVE} = $listRef

# Each element in $listRef is stored as $sequenceId|$peptideStartPosition. 
# The method also initializes $self->{CurrentTestType} to $POSITIVE (so the first peptides tested will be positive) 
# and $self->{CurrentTestIndex} to be 0.
#
# RETURN NULL
#######################################################################################################################################
sub initializePeptidePool{

    my ($self) = @_;

    my $peptidePool;

    my $model = $self->getModel();
    my $positives = $model->getPeptides($POSITIVE);
    
    #Add positives to pool
    my @positivePeptideList;
    
    foreach my $positiveSequenceId (keys %$positives){
	my $pOnePositions = $positives->{$positiveSequenceId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    
	    my $peptideCode = $positiveSequenceId  . "|" . $pOnePosition;
	    push (@positivePeptideList, $peptideCode);
	}
    }
    $peptidePool->{$POSITIVE} = \@positivePeptideList;

    #Add negatives to pool
    my $negatives = $model->getPeptides($NEGATIVE);
    my @negativePeptideList;
    
    foreach my $negativeSequenceId (keys %$negatives){
	my $pOnePositions = $negatives->{$negativeSequenceId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    
	    my $peptideCode = $negativeSequenceId  . "|" . $pOnePosition;
	    push (@negativePeptideList, $peptideCode);
	}
    }
    $peptidePool->{$NEGATIVE} = \@negativePeptideList;

    #Initialize internal data
    my $positiveCount = $model->getPeptideCount($POSITIVE);
    my $negativeCount = $model->getPeptideCount($NEGATIVE);    

    $self->writeLog("LeaveOneOutBenchmarker: loaded $positiveCount positives, $negativeCount negatives in Peptide Pool");

    $self->{PeptidePool} = $peptidePool;
    $self->{CurrentTestType} = $POSITIVE;
    $self->{CurrentTestIndex} = 0;
    

    $self->{TestSetPositiveCount} = $positiveCount;
    $self->{TestSetNegativeCount} = $negativeCount;

    if (($positiveCount + $negativeCount > 0)){
	#make sure at least one peptide was loaded (otherwise would cause problems if this method was called before model loaded peptides)
	
	$self->setPeptidePoolLoaded();
    }
}

sub getIterationCount{

    my ($self) = @_;

    my $model = $self->getModel();
    
    my $positiveCount = $model->getPeptideCount($POSITIVE);
    my $negativeCount = $model->getPeptideCount($NEGATIVE);

    my $totalCount = $negativeCount + $positiveCount;
    
    return $totalCount;
}





#######################################################################################################################################
# createTestSet
# Creates a test set composed of a single peptide (the one 'left out'). The peptide to be added is tracked by $self->{CurrentTestType}
# and $self->{CurrentTestIndex}, which are updated when createTrainingSet() is called. The peptide is taken from $self->{PeptidePool}.
#
# Adds the peptide to $self->{TestSet}; the peptide is represented as a dictionary of the form 

# $testPeptide->{$sequenceId}->{$pOnePosition} = $classification
#######################################################################################################################################
sub createTestSet{

    my ($self) = @_;

    return 0 if ($self->pipelineHasErrors());

    unless ($self->isPeptidePoolLoaded()){
	$self->initializePeptidePool();
    }

    my $currentTestType = $self->{CurrentTestType};
    my $model = $self->getModel();

    if ($currentTestType eq "Done"){
	$model->writeError("loo_iteration_error", "createTestSet", "Have tested all peptides but BenchmarkerPipeline is still iterating\n"); #tested
	return 0;
    }

    #get next peptide to be tested
    my $currentTestIndex = $self->{CurrentTestIndex};
    my $currentTestPeptideList = $self->{PeptidePool}->{$currentTestType};
    my $testPeptideCode = $currentTestPeptideList->[$currentTestIndex];
    
    #make sure the peptide hasn't been tested already
    if ($self->peptideWasTested($testPeptideCode, $currentTestType)){

	$model->writeError("loo_iteration_error", "createTestSet", "Attempting to test peptide $testPeptideCode classification $currentTestType a second time");
	#sort of tested -- gets caught by different error (testing on peptide in training set) -- which would always happen along with this error
	return 0;
    }
    else {
	$self->addTestedPeptide($testPeptideCode, $currentTestType);
    }

    #add peptide to list
    my ($sequenceId, $pOnePosition) = split('\|', $testPeptideCode);
    my $testPeptides;
    $testPeptides->{$sequenceId}->{$pOnePosition} = $currentTestType;
    $self->writeLog("LeaveOneOutBenchmarker: Adding peptide index $currentTestIndex type $currentTestType to test set"); #this might clutter up logs too much
    $self->{TestSet} = $testPeptides;
}


#######################################################################################################################################
# createTrainingSet
# Creates a training set consisting of all peptides in the system except for the one in the test set.
# 
# Peptides are added to $self->{TrainingSet} which is a dictionary of the form
# $self->{TrainingSet}->{$sequenceId}->{$pOnePosition} = $classification
#
# Both positives and negatives are added to the training set. This method also increments the current peptide that will go in the test
# set in the next iteration.
#######################################################################################################################################
sub createTrainingSet{

    my ($self) = @_;

    return 0 if ($self->pipelineHasErrors());

    my $trainingSet;

    my $currentTestType = $self->{CurrentTestType};
    my $currentTestIndex = $self->{CurrentTestIndex};

    my $peptidePool = $self->{PeptidePool};
    my $positives = $peptidePool->{$POSITIVE};

    #Add positives to the training set; skip $self->{CurrentTestIndex}
    my $i = -1;
    foreach my $positive(@$positives){
	$i++; #need to increment before skipping $currentTestIndex
	next if ($i == $currentTestIndex && $currentTestType == $POSITIVE);  #skip the one peptide in test set
	my ($sequenceId, $pOnePosition) = split('\|', $positive);
	$trainingSet->{$sequenceId}->{$pOnePosition} = $POSITIVE;
    }

    #Add negatives to the training set; skip $self->{CurrentTestIndex}
    my $model = $self->getModel();
    if ($model->trainsOnNegatives()){
	my $negatives = $peptidePool->{$NEGATIVE};
	my $tempCount = scalar(@$negatives);
	
	$i = -1;
	foreach my $negative(@$negatives){
	    $i++;
	    next if ($i == $currentTestIndex && $currentTestType == $NEGATIVE);
	    my ($sequenceId, $pOnePosition) = split('\|', $negative);
	    $trainingSet->{$sequenceId}->{$pOnePosition} = $NEGATIVE;
	}
    }
    
    $self->{TrainingSet} = $trainingSet;

    #Increment indices to the next peptide 
    my $currentTestPeptideList = $self->{PeptidePool}->{$currentTestType};
    my $typeTotalCount = scalar (@$currentTestPeptideList);

    if ($currentTestIndex == ($typeTotalCount - 1)){  #we have reached the end of the current peptide testing list
	if ($currentTestType eq $POSITIVE){ #done positives, switch to negatives
	    $self->{CurrentTestType} = $NEGATIVE;
	    $self->{CurrentTestIndex} = 0;
	}
	else { #done negatives, mark as Done
	    $self->{CurrentTestType} = "Done";
	}
    }
    else {  #just increment
	$currentTestIndex++;
	$self->{CurrentTestIndex} = $currentTestIndex;
    }

    my ($testPeptideSeqId, $testPeptideStartPosition,  $testPeptideClassification) = $self->getTestPeptide();

    #Check to make sure our test peptide is not in the training set. -- This assumes that for each iteration, createTestSet is called before createTrainingSet
    if ($self->{TrainingSet}->{$testPeptideSeqId}->{$testPeptideStartPosition}){
	my $classification =    $self->{TrainingSet}->{$testPeptideSeqId}->{$testPeptideStartPosition};
	my $errorMsg = "Added current test peptide to training set. Peptide is sequence $testPeptideSeqId position $testPeptideStartPosition.";
	my $currentTestType = $self->{CurrentTestType};
	$errorMsg .= "Currently testing peptides of type $currentTestType; test set peptide classification is $testPeptideClassification";
	$model->writeError("loo_iteration_error", "createTrainingSet", $errorMsg); #tested
    }
}


#######################################################################################################################################
# processResults
# Reads the results of applying the training set model to the test set. Only one score should be returned from the single peptide 
# in the test set; this is added to $self->{AllScores}->{$resultSetId}->{score}
# $resultSetId is described in SvmModel; it is a space-delimited string of the form
#
# "modbaseSeqId startPosition classification"
#
# If this is the last iteration of the system, $self->processAllResults() is called to finish up.
#######################################################################################################################################
sub processResults {

    my ($self, $resultFh) = @_;

    #Check if errors; if any exist, will be handled by Pipeline
    my $errors = $self->pipelineHasErrors();
    if ($errors){
	return $errors;
    }

    #Read $self->{TestSetIdList}; set by superclass when testing model
    my $testSetIdList = $self->{TestSetIdList};
    my $model = $self->getModel();
    
    #get score for single peptide
    my $resultSet = $model->parseResults($testSetIdList);

    my $size = scalar(@$testSetIdList);
    if ($size > 1){
	$model->writeError("loo_iteration_error", "processResults", "Have more than one peptide in test set (found $size peptides)");
	return 1;
    }

    my @resultSetKeys = keys %$resultSet;
    my $resultSetId = $resultSetKeys[0];
    my $score = $resultSet->{$resultSetId}->{score};

    $self->{AllScores}->{$resultSetId}->{score} = $score;

    if ($self->{CurrentTestType} eq "Done"){
	$self->processAllResults($resultFh);
    }

    return 0;
}

sub processAllResults{

    my ($self, $resultFh) = @_;

    my $model = $self->getModel();

    my $resultSet = $self->{AllScores};

    my $fpCount = 0;
    my $tpCount = 0;
    my $previousTpRate = 0.0;
    my $previousFpRate = 0.0;
    my $previousEvalue = -1;

    my $foundCriticalEvalue = 0;
    my $criticalTpRate = 0;
    my $criticalFpRate = 0;

    my $currentRunBins;

    my $truePositiveCount = $self->{TestSetPositiveCount};
    my $trueNegativeCount = $self->{TestSetNegativeCount};

    
    my $fpsAtEachTp = $self->{FpsAtEachTp};
    
    my $writeResultSet = 1;

    print $resultFh "0\t0\n" if ($writeResultSet == 1);
    foreach my $resultId (sort ({$resultSet->{$b}->{score} <=> $resultSet->{$a}->{score}} keys %$resultSet)){
		    
	if ($resultId =~ /$POSITIVE/){   
	    $tpCount++;
	    $fpsAtEachTp->{$tpCount}->{$fpCount}++;
	}
	if ($resultId =~ /$NEGATIVE/){
	    $fpCount++;
	}

	my $tpRate = ($tpCount * 1.0) / ($truePositiveCount * 1.0);
	my $fpRate = ($fpCount * 1.0) / ($trueNegativeCount * 1.0);

	my $eValue = $resultSet->{$resultId}->{score};

#	$currentRunBins = $self->binScore($tpCount, $fpCount, $tpRate, $fpRate, $eValue, $currentRunBins);
	print $resultFh "$fpRate\t$tpRate\t$eValue\t$resultId\t$fpCount\t$tpCount\n" if ($writeResultSet == 1);
	#change back if not work: all {score} to {eValue}; sort a before b
	if (($tpRate >= (1.0 - $fpRate)) && ($foundCriticalEvalue == 0)){
	    
	    #my ($criticalEvalue, $fpRateAtEvalue, $tpRateAtEvalue) = $self->interpolateEvalue($tpRate, $fpRate, $previousTpRate, $previousFpRate, $eValue, $previousEvalue);
	    $foundCriticalEvalue = 1;	    
	    $criticalTpRate = $tpRate;
	    $criticalFpRate = $fpRate;
	}
	
	$previousFpRate = $fpRate;
	$previousTpRate = $tpRate;
	$previousEvalue = $eValue;
    }
    print $resultFh "1\t1\n\n" if ($writeResultSet == 1);
    print $resultFh "Critical TP Rate: $criticalTpRate   FP: $criticalFpRate\n";
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



sub addTestedPeptide{
    my ($self, $testPeptideCode, $testPeptideType) = @_;
    $self->{TestedPeptides}->{$testPeptideType}->{$testPeptideCode} = 1;
}

sub peptideWasTested{
    my ($self, $testPeptideCode, $testPeptideType) = @_;
    if ($self->{TestedPeptides}->{$testPeptideType}->{$testPeptideCode} == 1){
	return 1;
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


sub getTestPeptide{
    my ($self) = @_;
    my $testPeptide = $self->{TestSet};
    my @testPeptideSeqIdList = keys %$testPeptide;
    my $testPeptideSeqId = $testPeptideSeqIdList[0];

    my $startPositions = $testPeptide->{$testPeptideSeqId};
    my @startPositionList = keys %$startPositions;
    my $startPosition = $startPositionList[0];
    
    my $classification = $testPeptide->{$testPeptideSeqId}->{$startPosition};
    return ($testPeptideSeqId, $startPosition, $classification);
}


sub isPeptidePoolLoaded{
    my ($self) = @_;
    return $self->{PeptidePoolLoaded};
}

sub setPeptidePoolLoaded{
    my ($self) = @_;
    $self->{PeptidePoolLoaded} = 1;
}



return 1;
