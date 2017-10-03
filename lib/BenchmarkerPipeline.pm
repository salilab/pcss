package BenchmarkerPipeline;

use strict;
use FileHandle;

use CleavageSiteModel;
use Benchmarker;
use JackknifeBenchmarker;
use LeaveOneOutBenchmarker;
use FastaReader;

use POSIX qw(ceil floor);
use List::Util qw(sum);



my $POSITIVE = "positive";
my $NEGATIVE = "negative";
my $TEST = "Test";
my $PROTEOME = "Proteome";


our @ISA = qw(Pipeline);


sub new{

    my ($class, $parameterFileName) = @_;

    my $self = $class->SUPER::new($parameterFileName, "model");

    bless ($self, $class);

    $self->initialize();
    $self->initializeBenchmarker();
    return $self;
}


sub initialize{
    my ($self) = @_;
    my $modelClass = $self->getParam("model_class");

    #create and initialize model object
    my $model = $modelClass->new();
    $self->setModel($model);
    $model->setParams($self->getParams());
    $model->setLog($self->{Log});
    $model->setResultFh($self->{ResultFh});
    $model->initialize();  
}



############################################################################################################################################################     
# loadPeptides()
# 
# Prepare all necesary input for loading peptides and pass to model for execution. Does both the loading and feature creation.
# Peptide information and features are specific to the type of model, so the model holds onto them after processing
# Return: NULL
#############################################################################################################################################################     
sub loadPeptides{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Loading peptides to prepare for benchmarking");

    #Get location of peptide pipeline run
    my $submitNodeSequencesDir = $self->getParam("head_node_preprocess_directory");   
    my $seqBatchTopDir = $self->getParam("seq_batch_top_directory");
    my $peptidePipelineResultsDir = $submitNodeSequencesDir . "/" . $seqBatchTopDir;    

    #Go through each seq_batch directory and get its peptide pipeline results file to read peptides
    my $seqBatchCount = $self->getParam("seq_batch_count");
    my $seqBatchDirectoryPrefix = $self->getParam("seq_batch_directory_prefix");
    my $model = $self->getModel();
    my $peptidePipelineResultFileName = $self->getParam("peptide_pipeline_result_file_name");

    my ($allSequenceCounts, $allPeptideCounts);

    for (my $i = 1; $i <= $seqBatchCount; $i++){
	my $fullSeqBatchDirectory = $peptidePipelineResultsDir . "/" . $seqBatchDirectoryPrefix . "_" . $i;
	my $fullPeptideResultsFile = $fullSeqBatchDirectory . "/" . $peptidePipelineResultFileName;
	
	#load peptides
	my @counts = $model->loadPeptides($fullPeptideResultsFile);

	#compile statistics
	if (@counts){
	    my $sequenceCounts = $counts[0];
	    my $peptideCounts = $counts[1];
	    foreach my $peptideType (keys %$peptideCounts){
		my $peptideCount = $peptideCounts->{$peptideType};
		my $sequenceCount = scalar(keys %{$sequenceCounts->{$peptideType}});
		$allPeptideCounts->{$peptideType} += $peptideCount;
		$allSequenceCounts->{$peptideType} += $sequenceCount;
	    }
	}
    }

    #Create SVM-readable features for all loaded peptides
    $model->createPeptideFeatures($POSITIVE);
    $model->createPeptideFeatures($NEGATIVE);

    foreach my $peptideType (keys %$allPeptideCounts){
	my $peptideCount = $allPeptideCounts->{$peptideType};
	my $sequenceCount = $allSequenceCounts->{$peptideType};
	$self->writeLog("BenchmarkerPipeline: Loaded $peptideCount peptides from $sequenceCount sequences with classification $peptideType");
    }

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "loadPeptides");
}


############################################################################################################################################################ 
# execute()
# 
# Run pipeline; creates the test and training sets, then trains and tests the model. In the webservice, only runs for one iteration (one iteration on each node).
# RETURN NULL
############################################################################################################################################################       
sub execute{
    
    my ($self) = @_;
    
    $self->loadPeptides();
    
    my $benchmarker = $self->getBenchmarker();
    
    $benchmarker->initialize();

    my $iterationCount = $benchmarker->getIterationCount();
    
    my $errors;
    
    $self->writeLog("BenchmarkerPipeline: running $iterationCount iterations");

    for (my $i = 0; $i < $iterationCount; $i++){
	if ($iterationCount % 100 == 0){
	    $self->writeLog("BenchmarkerPipeline: on interation $i");
	}
	$benchmarker->createTestSet();
	
	$benchmarker->createTrainingSet();
	
	$benchmarker->trainModel();
	
	$benchmarker->testModel();
	
	$errors = $benchmarker->processResults($self->{ResultFh});
	last if $errors;
    }    
    if ($errors){
	my $columnHeaderString = "TP_count\tFP_count\tScore";

	$self->writeNoOutput($errors, $columnHeaderString);
    }

}

############################################################################################################################################################     
# writeNoOutput
# 
# The input to this pipeline is the PeptidePipeline results file generated in the previous step. If there is no input due to no peptides being parsed from
# sequences in that step, or due to a fatal error in that step or in running this ApplicationPipeline, then indicate this to the results file here. In the case
# that the error was from the PeptidePipeline, it will have output an error message or keyword; grab that and write to the results file.
#
# PARAM  $noOutputReason: the keyword indicating no peptides were parsed or that there was an error
# PARAM  $columnHeaderString: tab delimited list of column names as read from $self->{ColumnInfo}
# RETURN NULL
############################################################################################################################################################     
sub writeNoOutput{
    my ($self, $noOutputReason, $columnHeaderString) = @_;

    $self->writeLog("BenchmarkerPipeline: No results will be output (reason: $noOutputReason). Writing results file indicating this.");

    my $resultFh = $self->{ResultFh};

    $self->writeResultFhGlobalError($noOutputReason, $columnHeaderString);
    
}




sub printPeptides{

    my ($self) = @_;
    my $model = $self->getModel();
    print STDERR "ALL POSITIVES:\n";
    $model->printPeptides($POSITIVE);
    print STDERR "ALL NEGATIVES:\n";
    $model->printPeptides($NEGATIVE);
}


####### Helper Methods #######

sub initializeBenchmarker{
    my ($self) = @_;

    #create and initialize benchmarker object
    my $benchmarkClass = $self->getParam("benchmark_class");  
    $benchmarkClass .= "Benchmarker";
    my $benchmarker = $benchmarkClass->new();
    $self->setBenchmarker($benchmarker);

    $benchmarker->setModel($self->{Model});
    $benchmarker->setLog($self->{Log});
    $benchmarker->setResultFh($self->{ResultFh});


}

sub setBenchmarker{
    my ($self, $benchmarker) = @_;
    $self->{Benchmarker} = $benchmarker;
}


sub getBenchmarker{
    my ($self) = @_;
    my $benchmarker = $self->{Benchmarker};
    return $benchmarker;
}

sub handleReaderError{
    my ($self, $inputFile) = @_;
    my $inputFh = FileHandle->new("<" . $inputFile);
    unless ($inputFh){
	print STDERR "BenchmarkerPipeline.pl: Exiting due to FastaReader error; could not open fasta file $inputFile for reading: $!\n";
	exit(1);
    }
    else {
	print STDERR "BenchmarkerPipeline.pl: Exiting due to FastaReader error; Reader did not return true value when reading $inputFile\n";
	#this should never happen, extreme sanity check
	exit(1);
    }
    
}



return 1;
