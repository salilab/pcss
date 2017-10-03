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

sub createPeptideIdentifierList{
    my ($self, $inputProteins, $exclude) = @_;

    my @peptideList;

    foreach my $seqId (keys %$inputProteins){
	my $startPositions = $inputProteins->{$seqId};
	foreach my $startPosition (keys %$startPositions){

	    #check exclude
	    next if ($exclude && $exclude->{$seqId}->{$startPosition});
			
	    #create unique identifier to be put into list for random selection
	    my $uniqueIdString = "$seqId" . "_" . $startPosition;
	    push (@peptideList, $uniqueIdString);
	}
    }
    return @peptideList;
}


#inputSeqs->{$sequenceId}->{$startPosition} = 1
sub createRandomPeptideList{

    my ($self, $inputSeqs, $count, $exclude) = @_;

    #convert peptides from my data structure to an array of strings uniquely specifying each one; throw out those specified in $exclude
    my @peptideList = $self->createPeptideIdentifierList($inputSeqs, $exclude);

    my $inputSize = scalar(@peptideList);
    
    my $randomSeqs;
    my $sampledEntries;
    my $sampledEntryCount = 1;
    my $sampleCount = 0;
    while ($sampledEntryCount <= $count){  

	#check count to avoid infinite loop
	$sampleCount++;
	if ($sampleCount > 1000000){
	    print STDERR "bailing out of random peptide list because sampled too much, infinite loop (benchmarker.pm)\n";   
	    #this avoids an infinite loop on the cluster if there are not enough peptides left to sample because of internal error.  
	    #TODO -- make nice log message that can easily be found in a search.  Considered sending email but that would probably be too complicated -- no place to write errors in output file in training mode
	    last;
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

sub trainModel{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Training model");

    my $model = $self->getModel();

    #create training set file name
    my $pipelineDir = $model->getParam("pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $trainingSetFileName = $model->getParam("svm_training_set_file_name");
    my $fullTrainingSetFileName = $pipelineDir . "/" . $runName . "/" . $trainingSetFileName;

    my $trainingSet = $self->getTrainingSet();

    #create full model file name
    my $pipelineDir = $model->getParam("pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $modelFileName = $model->getParam("svm_training_model_name");
    my $fullModelFile = $pipelineDir . "/" . $runName . "/" . $modelFileName;

    $model->trainModel($trainingSet, $fullTrainingSetFileName, $fullModelFile);

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "trainModel");

}

sub testModel{

    my ($self) = @_;
    my $model = $self->getModel();

    my $startTime = time();
    $self->writeLog("Testing model");
    

    #create test set file name
    my $pipelineDir = $model->getParam("pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $testSetFileName = $model->getParam("svm_test_set_file_name");
    my $fullTestSetFileName = $pipelineDir . "/" . $runName . "/" . $testSetFileName;

    my $testSet = $self->{TestSet};

    #create full model file name
    my $pipelineDir = $model->getParam("pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $modelFileName = $model->getParam("svm_training_model_name");
    my $fullModelFile = $pipelineDir . "/" . $runName . "/" . $modelFileName;
    
    my $testSetIdList = $model->testModel($testSet, $fullTestSetFileName, $fullModelFile);
    $self->{TestSetIdList} = $testSetIdList;

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "trainModel");

}


sub addToPeptideSet{

    my ($self, $finalPeptideSet, $peptideSet, $type) = @_;

    foreach my $peptideId (keys %$peptideSet){
	my $startPositions = $peptideSet->{$peptideId};
	
	foreach my $startPosition (keys %$startPositions){
	    $finalPeptideSet->{$peptideId}->{$startPosition} = $type;
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
