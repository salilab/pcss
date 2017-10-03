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

    my $self = $class->SUPER::new($parameterFileName);

    bless ($self, $class);

    return $self;
}



sub execute{

    my ($self) = @_;

    $self->loadPeptides();

    $self->createModel();

}



sub loadPeptides{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Loading peptides to prepare for benchmarking");

    #Get location of peptide pipeline run
    my $submitNodeSequencesDir = $self->getParam("head_node_preprocess_directory");   
    my $topLevelSeqDir = $self->getParam("top_level_seq_directory");
    my $peptidePipelineResultsDir = $submitNodeSequencesDir . "/" . $topLevelSeqDir;    

    #TODO -- change this so it reads the job name from the parameter file and uses it to go into preprocess (where results of peptide pipeline run are stored)
    #        right now the job name isn't being set until after the frontend calls submit.  

    #Get original input file (used in peptide pipeline) to get list of modbase seq ids (ignore protein sequence)
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $runName = $self->getParam("run_name");
    my $fastaFileName = $self->getParam("input_fasta_file_name");   
    my $fullFastaFileName = $pipelineDir . "/" . $runName . "/" . $fastaFileName;
    my $reader = FastaReader->new($fullFastaFileName, 1);
    $self->handleReaderError($fullFastaFileName) unless $reader;
    $reader->read();
    my $sequences = $reader->getSequences();

    my $model = $self->getModel();
    my $seqNotFoundKeyword = $self->getParam("keyword_no_modbase_for_uniprot");    
    my $peptidePipelineResultFileName = $self->getParam("peptide_pipeline_result_file_name");

    foreach my $sequenceInfo (keys %$sequences){

	#parse sequence header to get modbase seq id
	my @cols = split('\|', $sequenceInfo);
	my $modbaseSeqId = $cols[0];
	next if ($modbaseSeqId eq $seqNotFoundKeyword);
	
	#get peptide pipeline results for this sequence
	my $outputDir = $self->makeThreeLetterOutputDir($peptidePipelineResultsDir, $modbaseSeqId);
	my $peptideResultsFile = $outputDir . "/" . $peptidePipelineResultFileName;
	my $logMessage = $model->loadPeptides($peptideResultsFile);
    }
    $model->createPeptideFeatures($POSITIVE);
    $model->createPeptideFeatures($NEGATIVE);

    #take all peptides and put into training set
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
    my $pipelineDir = $model->getParam("pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $trainingSetFileName = $model->getParam("svm_complete_set_file_name");
    my $fullTrainingSetFileName = $pipelineDir . "/" . $runName . "/" . $trainingSetFileName;

    my $trainingSet = $self->{TrainingSet};

    #create full model file name
    my $pipelineDir = $model->getParam("pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $modelFileName = $model->getParam("svm_complete_model_name");
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
