package Benchmarker;

use JackknifeBenchmarker;
use LeaveOneOutBenchmarker;
use POSIX qw(ceil floor);
use strict;
use FileHandle;

sub new{
    my ($class, $model) = @_;
    my $self = {};
    bless $self, $class;
    $self->setModel($model);
    return $self;
}

#######################################################################################################################################
# trainModel()
# Uses model object to create an actual model by training on training-set peptides. Model will be applied to test set in next step.
#
# RETURN NULL
#######################################################################################################################################
sub trainModel{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Training model");

    return 0 if ($self->pipelineHasErrors());

    my $model = $self->getModel();

    #create training set file name
    my $pipelineDir = $model->getParam("cluster_pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $trainingSetFileName = $model->getParam("svm_training_set_file_name");
    my $fullTrainingSetFileName = $pipelineDir . "/" . $runName . "/" . $trainingSetFileName;

    my $trainingSet = $self->getTrainingSet();

    #create full model file name
    my $pipelineDir = $model->getParam("cluster_pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $modelFileName = $model->getParam("svm_training_model_name");
    my $fullModelFile = $pipelineDir . "/" . $runName . "/" . $modelFileName;

    #train
    $model->trainModel($trainingSet, $fullTrainingSetFileName, $fullModelFile);

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "trainModel");

}

#######################################################################################################################################
# testModel()
# Apply model generated in trainModel() to the test set, using model object.
# 
# RETURN NULL
#######################################################################################################################################
sub testModel{

    my ($self) = @_;
    my $model = $self->getModel();

    my $startTime = time();
    $self->writeLog("Testing model");
    
    return 0 if ($self->pipelineHasErrors());

    #create test set file name
    my $pipelineDir = $model->getParam("cluster_pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $testSetFileName = $model->getParam("svm_test_set_file_name");
    my $fullTestSetFileName = $pipelineDir . "/" . $runName . "/" . $testSetFileName;

    my $testSet = $self->{TestSet};

    #create full model file name
    my $pipelineDir = $model->getParam("cluster_pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $modelFileName = $model->getParam("svm_training_model_name");
    my $fullModelFile = $pipelineDir . "/" . $runName . "/" . $modelFileName;
    
    #testSetIdList is an array ref where each entry is of the format "$modbaseSequenceId $startPosition $status". 
    my $testSetIdList = $model->testModel($testSet, $fullTestSetFileName, $fullModelFile);
    $self->{TestSetIdList} = $testSetIdList;

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "trainModel");

}


#######################################################################################################################################
# createPeptideIdentifierList()
# Simple method that goes through a complex structure of peptides and returns an array where each entry is a unique id for a peptide.
# 
# PARAM  $inputProteins: peptides to examine (of the form $inputProteins->{$modbaseSeqId}->{$startPosition} = $seqInfo; $seqInfo unused here)
# PARAM  $exclude: peptides to exclude from the list, if any (same form as $inputProteins)
# RETURN @peptideList: array where each entry is of the form $modbaseSeqId_$startPosition
#######################################################################################################################################
sub createPeptideIdentifierList{
    my ($self, $inputProteins, $exclude) = @_;

    my @peptideList;

    foreach my $seqId (keys %$inputProteins){
	my $startPositions = $inputProteins->{$seqId};
	foreach my $startPosition (keys %$startPositions){

	    #check if we should exclude this peptide
	    next if ($exclude && $exclude->{$seqId}->{$startPosition});
			
	    #create unique identifier to be put into list for random selection
	    my $uniqueIdString = $seqId . "_" . $startPosition;
	    push (@peptideList, $uniqueIdString);
	}
    }
    return @peptideList;
}

#######################################################################################################################################
# createRandomPeptideList()
# Given a list of peptides, randomly select a subset of them. Note that whatever calls this method needs to make sure there are a number
# of peptides in $inputSeqs - $exclude that is greater than $count, to avoid inifite sampling (method checks this).
#
# PARAM  $inputSeqs: peptide pool from which to choose (of the form $inputSeqs->{$modbaseSeqId}->{$startPosition} = $seqInfo; $seqInfo unused)
# PARAM  $count: number of peptides to randomly choose
# PARAM  $exclude: peptides in here will be thrown back if selected. 
# RETURN $randomSeqs; hash of peptides that were selected (of the form $randomSeqs->{$modbaseSeqId}->{$startPosition} = 1
#######################################################################################################################################
sub createRandomPeptideList{

    my ($self, $inputSeqs, $count, $exclude) = @_;

    #put peptides in array to facilitate random selection
    my @peptideList = $self->createPeptideIdentifierList($inputSeqs, $exclude);

    my $inputSize = scalar(@peptideList);
    
    my $randomSeqs;
    my $sampledEntries;
    my $sampledEntryCount = 1;
    my $sampleCount = 0;
    $self->writeLog("Benchmarker: count = $count");
    while ($sampledEntryCount <= $count){  


	#check count to avoid infinite loop
	$sampleCount++;
	if ($sampleCount > 10000000){
	    $self->getModel()->writeError("infinite_sampling_loop", "createRandomPeptideList", "Bailing out of random peptide list due to infinite loop (sampled $sampleCount times)");
	    return 0;
	}

	#sample next entry without replacement
	my $nextSampledEntryIndex = int(rand($inputSize - 1));
	next if ($sampledEntries->{$nextSampledEntryIndex});   #if we already have this one, start over
	$sampledEntries->{$nextSampledEntryIndex} = 1;  
	
	#find peptide identifier at this index, add to list
	my $randomEntry = $peptideList[$nextSampledEntryIndex];
	my ($sequenceId, $startPosition) = split('\_', $randomEntry);
	$randomSeqs->{$sequenceId}->{$startPosition} = 1;  
	$sampledEntryCount++;
	
    }
    return $randomSeqs;
}



#######################################################################################################################################
# addToPeptideSet()
# Merge two sets of peptides and set the classification. 
# 
# PARAM  $finalPeptideSet: hash which may or may not be populated, of the form $finalPeptideSet->{$modbaseSeqId}->{$startPosition} = $type
# PARAM  $peptideSet: peptides in this set are added to $finalPeptideSet, of the from $finalPeptideSet->{$modbaseSeqId}->{$startPosition} = 1
# PARAM  $type: type that gets assigned to $finalPeptideSet
# RETURN $finalPeptideSet: same hash object as was input, in same format
#######################################################################################################################################
sub addToPeptideSet{

    my ($self, $finalPeptideSet, $peptideSet, $type) = @_;

    foreach my $modbaseSeqId (keys %$peptideSet){
	my $startPositions = $peptideSet->{$modbaseSeqId};
	
	foreach my $startPosition (keys %$startPositions){
	    $finalPeptideSet->{$modbaseSeqId}->{$startPosition} = $type;
	}
    }
    return $finalPeptideSet;
}



sub setLog{
    my ($self, $log) = @_;
    $self->{Log} = $log;
}

sub writeLog {

    my ($self, $message) = @_;
    my $logFh = $self->{Log};
    my ($sec,$min,$hour,$mday,$mon,$year,$wday, $yday,$isdst)=localtime(time);

    my $dateLine = sprintf ("%4d-%02d-%02d %02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

    print $logFh "$dateLine (Model):\t$message\n";
}


sub setResultFh{
    my ($self, $resultFh) = @_;
    $self->{ResultFh} = $resultFh;
}


sub setModel{
    my ($self, $model) = @_;
    $self->{Model} = $model;
}

sub getModel{
    my ($self) = @_;
    my $model = $self->{Model};
    return $model;
}


sub getTrainingSet{
    my ($self) = @_;
    my $trainingSet = $self->{TrainingSet};
    return $trainingSet;

}


sub getIterationCount{

    my ($self) = @_;
    print STDERR "method getIterationCount() should only be called from subclass of Benchmarker\n";
    exit(1);
}

sub createTestSet{

    my ($self) = @_;
    print STDERR "method createTestSet() should only be called from subclass of Benchmarker\n";
    exit(1);
}

sub createTrainingSet{

    my ($self) = @_;
    print STDERR "method createTrainingSet() should only be called from subclass of Benchmarker\n";
    exit(1);
}


sub processResults {

    my ($self, $resultFh) = @_;
    print STDERR "method processResults() should only be called from subclass of Benchmarker\n";
    exit(1);

}


sub processAllResults {

    my ($self, $resultFh) = @_;
    print STDERR "method processAllResults() should only be called from subclass of Benchmarker\n";
    exit(1);

}


sub initialize{
    my ($self) = @_;
    print STDERR "Method initialize() should only be called from subclass of Benchmarker.\n";
    exit(1);

}


sub writeDuration{
    my ($self, $startTime, $endTime, $functionName) = @_;

    my $difference = $endTime - $startTime;
    my $minutes = floor (($difference * 1.0 )/ 60.0);
    my $seconds = ($difference * 1.0) % 60.0;

    $self->writeLog("Duration of function $functionName: $minutes minutes $seconds seconds");
    if ($minutes > 60){
	$self->writeLog("Function $functionName took more than one hour to run");
    }
}



return 1;
