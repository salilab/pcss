package ApplicationPipeline;

use strict;
use FileHandle;

use CleavageSiteModel;
use Pipeline;

use POSIX qw(ceil floor);
use List::Util qw(sum);

my $POSITIVE = "positive";
my $NEGATIVE = "negative";
my $TEST = "Test";
my $PROTEOME = "Application";


our @ISA = qw(Pipeline);

####### Top level methods (called by pipeline script) #######
sub new{

    my ($class, $parameterFileName) = @_;

    my $self = $class->SUPER::new($parameterFileName);

    bless ($self, $class);

    return $self;
}

sub loadPeptides{

    my ($self) = @_;
    my $model = $self->getModel();
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $runName = $self->getParam("run_name");
    my $peptideFile = $self->getParam("peptide_pipeline_result_file_name");  
    
    my $fullPeptideFileName = $pipelineDir . "/" . $runName . "/" . $peptideFile;

    my $logMessage = $model->loadPeptides($fullPeptideFileName, $PROTEOME);
    
    $model->createPeptideFeatures($PROTEOME);  #TODO - this type should probably be parameterized since it is set initially elsewhere in the server

    $self->writeLog($logMessage);
}

sub execute{

    my ($self) = @_;

    $self->loadColumnInfo();

    $self->loadPeptides();

    if ($self->checkInput() == 0){

	$self->createApplicationSet();
	
	$self->applyModel();
	
	$self->processResults();
    }
    else {
	$self->writeLog("ApplicationPipeline: model reports no peptides parsed");
	$self->writeNoInput();
    }
}


####### Pipeline secondary level methods  #######

sub createApplicationSet{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Creating Application Set");

    my $model = $self->getModel();

    my $peptides = $model->getPeptides($PROTEOME);
    
    my $applicationSet;

    foreach my $peptideSequenceId (keys %$peptides){
	my $pOnePositions = $peptides->{$peptideSequenceId};
	
	foreach my $pOnePosition (keys %$pOnePositions){

	    $applicationSet->{$peptideSequenceId}->{$pOnePosition} = $PROTEOME;
	}
    }
  
    $self->{ApplicationSet} = $applicationSet;
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "createApplicationSet");
}


sub applyModel{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Applying model");
    
    my $model = $self->getModel();
    my $applicationSet = $self->{ApplicationSet};
    my $pipelineDir = $model->getParam("pipeline_directory");
    my $runName = $model->getParam("run_name");

    my $applicationSetFileName = $self->getParam("svm_application_set_file_name");
    my $fullApplicationSetFileName = $pipelineDir . "/" . $runName . "/" . $applicationSetFileName;

    my $modelFileName = $self->getParam("svm_application_model");

    my $applicationSetIdList = $model->testModel($applicationSet, $fullApplicationSetFileName, $modelFileName);
    $self->{ApplicationSetIdList} = $applicationSetIdList;

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "applyModel");

}



sub checkInput{
    my ($self) = @_;
    my $model = $self->getModel();
    if ($model->getNoPeptidesParsed() == 1){
	return 1;
    }
    else {
	return 0;
    }

}

sub writeNoInput{
    my ($self) = @_;

    #get original input file
    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $pipelineInputFile = $self->getParam("peptide_pipeline_result_file_name");
    my $fullPipelineInputFile = $pipelineDir . "/" . $run . "/" . $pipelineInputFile;

    my $pipelineInputFh = FileHandle->new("<" . $fullPipelineInputFile) || die "could not open $fullPipelineInputFile\n";
    
    my $resultFh = $self->{ResultFh};

    #get column names
    my $line = <$pipelineInputFh>;
    chomp $line;
    print $resultFh "SVM Score\t$line\n";

    #write no result keyword
    my $noPeptidesParsedKeyword = $self->getParam("keyword_no_peptides_parsed");
    print $resultFh $noPeptidesParsedKeyword . "\n";
}


sub processResults {

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Processing model application results");


    my $model = $self->getModel();

    my $applicationSetIdList = $self->{ApplicationSetIdList};
    my $resultSet = $model->parseResults($applicationSetIdList);
    $model->addResultsToPeptideInfo($resultSet, $PROTEOME);
    $model->addBenchmarkScoresToPeptideInfo($PROTEOME);

    #add svm results to original pipeline result file (that itself was used as input to the SVM)
    my $resultCount = scalar (keys %$resultSet);

    #get original input file
    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $pipelineInputFile = $self->getParam("peptide_pipeline_result_file_name");
    my $fullPipelineInputFile = $pipelineDir . "/" . $run . "/" . $pipelineInputFile;

    my $pipelineInputFh = FileHandle->new("<" . $fullPipelineInputFile) || die "could not open $fullPipelineInputFile\n";
    
    my $resultFh = $self->{ResultFh};

    my $columnMap;
    my $lineCounter = 0;

    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");
    my $keywordApplicationDisplay = $self->getParam("keyword_application_display");

    my $columnHeaderString = $self->makeColumnHeaderString();
    print $resultFh $columnHeaderString . "\n";
    
    my $peptides = $model->getPeptides($PROTEOME);
    my $columns = $self->{ColumnInfo};

    foreach my $peptideSequenceId (keys %$peptides){
	my $pOnePositions = $peptides->{$peptideSequenceId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $outputLine = "";
	    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
		my $modes = $columns->{$columnShortName}->{modes};   
		if ($modes =~ /$keywordFeatureInternal/ || $modes =~ /$keywordApplicationDisplay/){
		    
		    my $columnName = $columns->{$columnShortName}->{displayName};
		    my $method = $columns->{$columnShortName}->{method};
		    if ($method eq "universal"){
		    
			my $nextOutputValue = $seqInfo->{$columnName};
			$outputLine .= $nextOutputValue . "\t";
		    }
		    else {
			if ($self->runMethod($method)){
			    my $nextOutputValue = $seqInfo->{$columnName};
			    $outputLine .= $nextOutputValue . "\t";
			}
		    }
		}
	    }	    
	    print $resultFh $outputLine . "\n"; 
	}
    }
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "processResults");
}

sub runMethod{

    my ($self, $methodName) = @_;
    
    my $runMethodParamValue = $self->getParam($methodName);

    if ($runMethodParamValue eq "yes"){
	return 1;
    }
    else {
	return 0;
    }
}


sub makeColumnHeaderString{

    my ($self) = @_;

    my $columnHeaderString = "";
    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");
    my $keywordApplicationDisplay = $self->getParam("keyword_application_display");
    
    my $columns = $self->{ColumnInfo};
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   
	if ($modes =~ /$keywordFeatureInternal/ || $modes =~ /$keywordApplicationDisplay/){
		    
	    my $displayName = $columns->{$columnShortName}->{displayName};
	    my $method = $columns->{$columnShortName}->{method};
	    if ($method eq "universal"){
		$columnHeaderString .= $displayName . "\t";
	    }
	    else {
		if ($self->runMethod($method)){
		    $columnHeaderString .= $displayName . "\t";
		}
	    }
	}
    }
    return $columnHeaderString;
}


sub loadColumnInfo{

    my ($self) = @_;
    my $columnInfoFile = $self->getParam("column_info_file");
    my $fh = FileHandle->new("<" . $columnInfoFile) || die "could not open column info file $columnInfoFile\n";

    my $columnInfo;
    my $counter = 1;
    while (<$fh>){
	chomp;
	my $line = $_;
	my ($columnName, $columnShortName, $columnModes, $columnMethod, $columnDescription) = split('\t', $line);

	$columnInfo->{$columnShortName}->{displayName} = $columnName;
	$columnInfo->{$columnShortName}->{modes} = $columnModes;
	$columnInfo->{$columnShortName}->{method} = $columnMethod;
	$columnInfo->{$columnShortName}->{description} = $columnDescription;
	$columnInfo->{$columnShortName}->{displayOrder} = $counter;
	$counter++;
    }
    $self->{ColumnInfo} =  $columnInfo;
}


return 1;
