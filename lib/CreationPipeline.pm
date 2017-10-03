package CreationPipeline;

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


####### Top level methods (called by pipeline script) #######
sub new{

    my ($class, $parameterFileName) = @_;

    my $self = $class->SUPER::new($parameterFileName, "model");

    bless ($self, $class);

    $self->initialize();

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




sub execute{

    my ($self) = @_;

    $self->loadPeptides();
    $self->createModel();
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

    my $trainingSet; 
    my $positives = $model->getPeptides($POSITIVE);
    $trainingSet = $self->addToPeptideSet($trainingSet, $positives, $POSITIVE);
    my $negatives = $model->getPeptides($NEGATIVE);
    $trainingSet = $self->addToPeptideSet($trainingSet, $negatives, $NEGATIVE);
    $self->{TrainingSet} = $trainingSet;

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "loadPeptides");
}


sub createModel{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Creating model");

    my $model = $self->getModel();

    #create training set file name
    my $pipelineDir = $model->getParam("cluster_pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $trainingSetFileName = $model->getParam("svm_complete_set_file_name");
    my $fullTrainingSetFileName = $pipelineDir . "/" . $runName . "/" . $trainingSetFileName;

    my $trainingSet = $self->{TrainingSet};

    #create full model file name
    my $pipelineDir = $model->getParam("cluster_pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $modelFileName = $model->getParam("user_created_svm_model_name");
    my $fullModelFile = $pipelineDir . "/" . $runName . "/" . $modelFileName;

    $model->trainModel($trainingSet, $fullTrainingSetFileName, $fullModelFile);

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "createModel");

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

sub makeThreeLetterOutputDir{
    my ($self, $sequencesDir, $modbaseSeqId) = @_;
    $modbaseSeqId =~ /^(\S\S\S).*/;
    my $threeLetterDir = $1;
    
    my $outputDir = $sequencesDir . "/" . $threeLetterDir . "/" . $modbaseSeqId;

    return $outputDir;
}
