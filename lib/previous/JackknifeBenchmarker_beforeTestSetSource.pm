package JackknifeBenchmarker;

use strict;
use FileHandle;



use POSIX qw(ceil floor);
use List::Util qw(sum);


my $POSITIVE = "positive";
my $NEGATIVE = "negative";
my $TEST = "Test";
my $PROTEOME = "Proteome";

our @ISA = qw(Benchmarker);

sub new {
    my ($class, $model) = @_;
    my $self = $class->SUPER::new($model);

    bless ($self, $class);
   
    return $self;
}

sub getIterationCount{

    my ($self) = @_;

    my $model = $self->getModel();
    my $iterationCount = $model->getParam("iteration_count");

    return $iterationCount;
}

sub createTestSet{

    my ($self) = @_;

    my $model = $self->getModel();

    my $positives = $model->getPeptides($POSITIVE);
    my $positiveCount = $model->getPeptideCount($POSITIVE);
    my $testSetPercentage = $model->getParam("test_set_percentage");
    my $testSetPositiveCount = floor(($positiveCount * 1.0) * ($testSetPercentage * 1.0));

    my $negativeTestSetRatio = $model->getParam("negative_test_set_ratio");
    my $testSetNegativeCount = $testSetPositiveCount * $negativeTestSetRatio;

    my $testSetSource = $model->getParam("test_set_source");

    print STDERR "creating test set:  total positive count $positiveCount percentage: $testSetPercentage test set positive count: $testSetPositiveCount negative count; $testSetNegativeCount source $testSetSource\n";

    my $testSetPositives = $self->createRandomPeptideList($positives, $testSetPositiveCount, 0, $testSetSource);

    
    my $negatives = $model->getPeptides("Negative");
    my $testSetNegatives = $self->createRandomPeptideList($negatives, $testSetNegativeCount, 0, $testSetSource);
    
    $self->{TestSetPositiveCount} = $testSetPositiveCount;
    $self->{TestSetNegativeCount} = $testSetNegativeCount;

    my $finalTestSet;
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetPositives, $POSITIVE);
    $finalTestSet = $self->addToPeptideSet($finalTestSet, $testSetNegatives, $NEGATIVE);


    $self->{TestSet} = $finalTestSet;
}

sub createTrainingSet{

    my ($self) = @_;

    my $sourceDict = $self->buildSourceHash("training_set_source");
    my $testSet = $self->{TestSet};    
    my $model = $self->getModel();

    
    my $negatives = $model->getPeptides($NEGATIVE);

    print STDERR "training set positives\n";
    #positives
    my $positives = $model->getPeptides($POSITIVE);
    my $source = $model->getParam("training_set_source");
    my $positiveTrainingSetCount = 0;
    my $trainingSet;
    my $initialPositiveCount = $model->getPeptideCount($POSITIVE);
    print STDERR "initial positive count: $initialPositiveCount\n";
    
    foreach my $positiveSequenceId (keys %$positives){
	my $pOnePositions = $positives->{$positiveSequenceId};
	
	foreach my $pOnePosition (keys %$pOnePositions){

	    if ($source && $source ne "none"){   #if source was selected, enforce that here
		my $peptideSource = $model->getPeptideSource($positiveSequenceId, $pOnePosition);
		next unless ($sourceDict->{$peptideSource});
	    }
	    unless ($testSet->{$positiveSequenceId}->{$pOnePosition}){   #check if this cleavage sequence is in the test set
		$trainingSet->{$positiveSequenceId}->{$pOnePosition} = $POSITIVE;
		$positiveTrainingSetCount++;
	    }
	}
    }
    print STDERR "added $positiveTrainingSetCount positives to training set\n";
    #negatives
    my $trainsOnNegatives = $model->trainsOnNegatives();
    if ($trainsOnNegatives == 1){
	
	my $negativeTrainingSetRatio = $model->getParam("negative_training_set_ratio");
	my $negativeTrainingSetCount = $negativeTrainingSetRatio * $positiveTrainingSetCount;
	print STDERR "adding $negativeTrainingSetCount negatives to training set\n";
	my $trainingSetNegatives = $self->createRandomPeptideList($negatives, $negativeTrainingSetCount, $testSet, $source);

	$trainingSet = $self->addToPeptideSet($trainingSet, $trainingSetNegatives, $NEGATIVE);

    }
    $self->{TrainingSet} = $trainingSet;
}

sub processAllResults{
    my ($self, $resultFh) = @_;
    print STDERR "jackknife: processing allr esults\n";
    my $fpsAtEachTp = $self->{FpsAtEachTp};
    my @criticalEvalues = @{$self->{CriticalEValues}};
    my @tpRatesAtEvalue = @{$self->{tpRatesAtEvalue}};
    my @fpRatesAtEvalue = @{$self->{fpRatesAtEvalue}};

    my $truePositiveCount = $self->{TestSetPositiveCount};
    my $trueNegativeCount = $self->{TestSetNegativeCount};
   
    my $resultSetCount = scalar(@criticalEvalues);
  
    print $resultFh "fpr\ttpr\n";
    print $resultFh "0.0\t0.0\n";
    my $totalFps = $trueNegativeCount * $resultSetCount;
    foreach my $tp (sort ({$a <=> $b} keys %$fpsAtEachTp)){
	my $fpCounts = $fpsAtEachTp->{$tp};
	my $totalFpsAtThisTp = 0;

	my @allFps;
	
	foreach my $fpCount (keys %$fpCounts){
	    my $fpCountValue = $fpCounts->{$fpCount};
	    $fpCount *= $fpCountValue;
	    $totalFpsAtThisTp += $fpCount;

	    my $fpRate = ($fpCount * 1.0)  / ($totalFps * 1.0);  
	    push (@allFps, $fpRate);  #prepare standard deviation list
	}
	my $fpRateAtThisTp = ($totalFpsAtThisTp * 1.0) / ($totalFps * 1.0);
	my $tpRate = ($tp * 1.0) / ($truePositiveCount * 1.0);


	#get standard deviation
	my $squaredFpSum = 0;
       
	foreach my $fpRate (@allFps){

	    my $squaredFp = $fpRate * $fpRate;
	    $squaredFpSum += $squaredFp;
	}
	my $count = scalar(@allFps);
	my $tempFpAverage = ($squaredFpSum * 1.0) / ($count * 1.0);
	$squaredFpSum /= $count;
	my $squaredAverage = $tempFpAverage * $tempFpAverage;
	my $stdDev = sqrt($squaredFpSum - $squaredAverage);

	print $resultFh "$fpRateAtThisTp\t$tpRate\t$stdDev\n";
    }
       
    print $resultFh "1.0\t1.0\n\n\n";

    print $resultFh "***False Positive Rates, True Positive Rates at each SVM score:\n";

#    $self->writeScoreBinRates($trueNegativeCount, $truePositiveCount, $resultSetCount, $resultFh);

    my $totalEvalue = 0;
    foreach my $eValue (@criticalEvalues){
	
	$totalEvalue += $eValue;
    }
    my $avgEvalue = $totalEvalue / ($resultSetCount * 1.0);
    print $resultFh "\naverage critical evalue: $avgEvalue\n";


    my $totalTpAtEvalue = 0;
    foreach my $tpRate (@tpRatesAtEvalue){
	$totalTpAtEvalue += $tpRate;
    }
    my $avgTpAtEvalue = $totalTpAtEvalue / ($resultSetCount * 1.0);
    print $resultFh "\naverage tp rate at critical eValue: $avgTpAtEvalue\n";

    my $totalFpAtEvalue = 0;
    foreach my $fpRate (@fpRatesAtEvalue){
	$totalFpAtEvalue += $fpRate;
    }
    my $avgFpAtEvalue = $totalFpAtEvalue / ($resultSetCount * 1.0);
    print $resultFh "\naverage fp rate at critical eValue: $avgFpAtEvalue\n";
}


sub processResults {

    my ($self, $resultFh) = @_;
    
    my $testSetIdList = $self->{TestSetIdList};
    my $model = $self->getModel();
    my $resultSet = $model->parseResults($testSetIdList);

    my $fpCount = 0;
    my $tpCount = 0;
    my $previousTpRate = 0.0;
    my $previousFpRate = 0.0;
    my $previousEvalue = -1;

    my $foundCriticalEvalue = 0;
    my $currentRunBins;

    my $truePositiveCount = $self->{TestSetPositiveCount};
    my $trueNegativeCount = $self->{TestSetNegativeCount};

    my $fpsAtEachTp = $self->{FpsAtEachTp};
    
    my $writeResultSet = 1;

    print $resultFh "Results for first iteration:\n" if ($writeResultSet == 1);
    print $resultFh "0\t0\n" if ($writeResultSet == 1);
    foreach my $resultId (sort ({$resultSet->{$b}->{score} <=> $resultSet->{$a}->{score}} keys %$resultSet)){
		    
	if ($resultId =~ /Positive/){   
	    $tpCount++;
	    $fpsAtEachTp->{$tpCount}->{$fpCount}++;
	}
	if ($resultId =~ /Negative/){
	    $fpCount++;
	}

	my $tpRate = ($tpCount * 1.0) / ($truePositiveCount * 1.0);
	my $fpRate = ($fpCount * 1.0) / ($trueNegativeCount * 1.0);
	my $eValue = $resultSet->{$resultId}->{score};
#	$currentRunBins = $self->binScore($tpCount, $fpCount, $tpRate, $fpRate, $eValue, $currentRunBins);
	print $resultFh "$fpRate\t$tpRate\t$eValue\t\t$resultId\t$fpCount\t$tpCount\n" if ($writeResultSet == 1);
	#change back if not work: all {score} to {eValue}; sort a before b
	if (($tpRate >= (1.0 - $fpRate)) && ($foundCriticalEvalue == 0)){
	    
	    #my ($criticalEvalue, $fpRateAtEvalue, $tpRateAtEvalue) = $self->interpolateEvalue($tpRate, $fpRate, $previousTpRate, $previousFpRate, $eValue, $previousEvalue);
	    $foundCriticalEvalue = 1;	    
	    
	    push (@{$self->{CriticalEValues}}, $eValue);
	    push (@{$self->{tpRatesAtEvalue}}, $tpRate);
	    push (@{$self->{fpRatesAtEvalue}}, $fpRate);
	}
	
	$previousFpRate = $fpRate;
	$previousTpRate = $tpRate;
	$previousEvalue = $eValue;
    }
    print $resultFh "1\t1\n\n" if ($writeResultSet == 1);
   
}

sub initialize{

    my ($self) = @_;
}



return 1;
