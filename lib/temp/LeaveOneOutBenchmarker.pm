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
   
    return $self;
}

sub initialize{
    my ($self) = @_;
   
    $self->initializePeptidePool(); 


}

#create peptide pool that represents the only peptides from which the training set and test set will be drawn
sub initializePeptidePool{

    my ($self) = @_;

    my $peptidePool;

    my $model = $self->getModel();
    my $positives = $model->getPeptides($POSITIVE);
    
    #positives
    my @positivePeptideList;
    
    foreach my $positiveSequenceId (keys %$positives){
	my $pOnePositions = $positives->{$positiveSequenceId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    
	    my $peptideCode = $positiveSequenceId  . "|" . $pOnePosition;
	    push (@positivePeptideList, $peptideCode);
	}
    }
    $peptidePool->{$POSITIVE} = \@positivePeptideList;

    #negatives
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


    my $positiveCount = $model->getPeptideCount($POSITIVE);
    my $negativeCount = $model->getpeptideCount($NEGATIVE);    

    $self->writeLog("leave one out benchmarker: have $positiveCount positives, $negativeCount negatives in training set");

    $self->{PeptidePool} = $peptidePool;
    $self->{CurrentTestType} = $POSITIVE;
    $self->{CurrentTestIndex} = 0;
    

    $self->{TestSetPositiveCount} = $positiveCount;
    $self->{TestSetNegativeCount} = $negativeCount;

}

sub getIterationCount{

    my ($self) = @_;

    my $model = $self->getModel();
    
    my $positiveCount = $model->getPeptideCount($POSITIVE);

    my $negativeTestSetRatio = $model->getParam("negative_test_set_ratio");
    
    my $negativeCount = $positiveCount * $negativeTestSetRatio;
    my $totalCount = $negativeCount + $positiveCount;
    
    return $totalCount;
}


sub createTestSet{

    my ($self) = @_;

    my $currentTestType = $self->{CurrentTestType};

    if ($currentTestType eq "Done"){
	print STDERR "error: ran out of test peptides but still iterating\n";
	exit(1);
    }

    #get next peptide to be tested
    my $currentTestIndex = $self->{CurrentTestIndex};
    my $currentTestPeptideList = $self->{PeptidePool}->{$currentTestType};

    my $testPeptideCode = $currentTestPeptideList->[$currentTestIndex];
    my ($sequenceId, $pOnePosition) = split('\|', $testPeptideCode);
    my $testPeptides;
    $testPeptides->{$sequenceId}->{$pOnePosition} = $currentTestType;
    $self->writeLog("adding peptide index $currentTestIndex type $currentTestType to test set");
    $self->{TestSet} = $testPeptides;
}


sub createTrainingSet{

    my ($self) = @_;

    my $trainingSet;

    my $currentTestType = $self->{CurrentTestType};
    my $currentTestIndex = $self->{CurrentTestIndex};

    my $peptidePool = $self->{PeptidePool};
    my $positives = $peptidePool->{$POSITIVE};

    my $i = -1;
    #positives
    foreach my $positive(@$positives){
	$i++;

	next if ($i == $currentTestIndex && $currentTestType == $POSITIVE);  #skip the one peptide in test set

	my ($sequenceId, $pOnePosition) = split('\|', $positive);
	$trainingSet->{$sequenceId}->{$pOnePosition} = $POSITIVE;
    }

    #negatives
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

    #increment indices to the next peptide after that
    my $currentTestPeptideList = $self->{PeptidePool}->{$currentTestType};
    my $typeTotalCount = scalar (@$currentTestPeptideList);
    if ($currentTestIndex == ($typeTotalCount - 1)){  #we have reached the end of the current peptide testing list
	if ($currentTestType eq $POSITIVE){
	    $self->{CurrentTestType} = $NEGATIVE;
	    $self->{CurrentTestIndex} = 0;
	}
	else {
	    $self->{CurrentTestType} = "Done";
	}
    }
    else {
	$currentTestIndex++;
	$self->{CurrentTestIndex} = $currentTestIndex;
    }
}



sub processResults {

    my ($self, $resultFh) = @_;

    my $testSetIdList = $self->{TestSetIdList};
    my $model = $self->getModel();
    my $resultSet = $model->parseResults($testSetIdList);

    my $size = scalar(@$testSetIdList);
    print STDERR "sanity check: test set id list size: $size\n";

    my @resultSetKeys = keys %$resultSet;
    my $resultSetId = $resultSetKeys[0];
    my $score = $resultSet->{$resultSetId}->{score};

    $self->{AllScores}->{$resultSetId}->{score} = $score;
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



return 1;
