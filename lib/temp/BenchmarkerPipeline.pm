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

####### Top level methods (called by pipeline script) #######
sub new{

    my ($class, $parameterFileName) = @_;

    my $self = $class->SUPER::new($parameterFileName);

    bless ($self, $class);

    $self->initializeBenchmarker();
    return $self;
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

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "loadPeptides");

}



sub execute{

    my ($self) = @_;

    $self->loadPeptides();

    my $benchmarker = $self->getBenchmarker();

#    my $iterationCount = $benchmarker->getIterationCount();

#    for (my $i = 0; $i < $iterationCount; $i++){

	$benchmarker->createTestSet();

	$benchmarker->createTrainingSet();
	
	$benchmarker->trainModel();
	
	$benchmarker->testModel();
	
	$benchmarker->processResults($self->{ResultFh});
#    }    
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

sub makeThreeLetterOutputDir{
    my ($self, $sequencesDir, $modbaseSeqId) = @_;
    $modbaseSeqId =~ /^(\S\S\S).*/;
    my $threeLetterDir = $1;
    
    my $outputDir = $sequencesDir . "/" . $threeLetterDir . "/" . $modbaseSeqId;

    return $outputDir;
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
