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

############################################################################################################################################################     
# execute
#
# Handle feature creation and model application. Checks to make sure no errors occurred previously, or if there was no input to process
# (writes appropriate output if so)
# 
# RETURN: NULL
############################################################################################################################################################     
sub execute{

    my ($self) = @_;

    $self->loadPeptides();
    	
    $self->createApplicationSet();
    
    $self->applyModel();
    
    $self->processResults();
    
    $self->writeResults();
}


############################################################################################################################################################     
# loadPeptides
# 
# Prepare all necesary input for loading peptides and pass to model for execution. Does both the loading and feature creation.
# Peptide information and features are specific to the type of model, so the model holds onto them after processing
# Return: NULL
#############################################################################################################################################################     
sub loadPeptides{

    my ($self) = @_;

    my $startTime = time();

    $self->writeLog("ApplicationPipeline: loading peptides");

    my $model = $self->getModel();

    #get full peptide results file name
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $runName = $self->getParam("run_name");
    my $peptideFile = $self->getParam("peptide_pipeline_result_file_name");  
    my $fullPeptideFileName = $pipelineDir . "/" . $runName . "/" . $peptideFile;

    #load peptides and create features
    my @counts = $model->loadPeptides($fullPeptideFileName, $PROTEOME);
    if (@counts){
	my $sequenceCounts = $counts[0];
	my $peptideCounts = $counts[1];
	foreach my $peptideType (keys %$peptideCounts){
	    my $peptideCount = $peptideCounts->{$peptideType};
	    my $sequenceCount = scalar(keys %{$sequenceCounts->{$peptideType}});
	    $self->writeLog("ApplicationPipeline: Loaded $peptideCount peptides from $sequenceCount sequences with classification $peptideType");
	}
    
	$model->createPeptideFeatures($PROTEOME);  #TODO - this type should probably be parameterized since it is set initially elsewhere in the server
    }
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "loadPeptides");
}


############################################################################################################################################################     
# createApplicationSet
# 
# Get all peptides that have been loaded by the model and put them in $self->{ApplicationSet}. Peptides are formed in a dictionary with the following
# format: $applicationSet->{$modbaseSequenceId}->{$peptideStartPosition} = $PROTEOME.
#
# $self->{ApplicationSet} is thus a quick way to get a handle on the list of peptides to be processed by this pipeline. This is slightly heavy-handed,
# since we are passing all these IDs around in different data structures; however, it is necessary to remain consistent with training/test mode (model application
# and model testing on a test set are the same procedure).
#
# RETURN: NULL
############################################################################################################################################################     
sub createApplicationSet{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Creating application set to be processed in this pipeline run");

    return 0 if ($self->pipelineHasErrors());

    my $model = $self->getModel();
    my $peptides = $model->getPeptides($PROTEOME);
    my $applicationSet;

    #iterate through all peptides in model of type $PROTEOME (should be everything loaded); add to $applicationSet
    foreach my $peptideSequenceId (keys %$peptides){
	my $startPositions = $peptides->{$peptideSequenceId};
	
	foreach my $startPosition (keys %$startPositions){

	    $applicationSet->{$peptideSequenceId}->{$startPosition} = $PROTEOME;
	}
    }
    $self->{ApplicationSet} = $applicationSet;
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "createApplicationSet");
}

############################################################################################################################################################     
# applyModel
#
# Prepare input which model will process and pass to model for execution. Model will write file in its own format for processing and then actually apply
# itself to that file.
#
# Model writes scoring results into a file without labeling them; the only way we know which score is for which peptide is by saving the order in which 
# they were written to the file, which is stored in $self->{ApplicationSetIdList} and matched later when results are processed
#
# RETURN: NULL
############################################################################################################################################################     
sub applyModel{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Applying model");
    
    return 0 if ($self->pipelineHasErrors());

    my $model = $self->getModel();
    my $applicationSet = $self->{ApplicationSet};

    #get name of file that model will write
    my $pipelineDir = $model->getParam("cluster_pipeline_directory");
    my $runName = $model->getParam("run_name");
    my $applicationSetFileName = $self->getParam("svm_application_set_file_name");
    my $fullApplicationSetFileName = $pipelineDir . "/" . $runName . "/" . $applicationSetFileName;

    #$modelFileName is the model which was trained and benchmarked by the server in a separate training run
    my $modelFileName = $self->getParam("svm_application_model");

    #applicationSetIdList is an array ref where each entry is of the format "$modbaseSequenceId $startPosition $status". 
    my $applicationSetIdList = $model->testModel($applicationSet, $fullApplicationSetFileName, $modelFileName); 
    $self->{ApplicationSetIdList} = $applicationSetIdList;

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "applyModel");

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


sub setErrors{
    my ($self, $errorKeyword) = @_;
    $self->{errors} = $errorKeyword;
}

sub getErrors{
    my ($self) = @_;
    my $errors = $self->{errors};
    return $errors;
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

    $self->writeLog("ApplicationPipeline: No results will be output (reason: $noOutputReason). Writing results file indicating this.");

    my $resultFh = $self->{ResultFh};

    $self->writeResultFhGlobalError($noOutputReason, $columnHeaderString);
}



############################################################################################################################################################     
# processResults
#
# Given the ApplicationSet loaded previously, coordinate with the model to assign scores generated in applyModel() to the peptides contained in the set.
# When done, write results to $self->{ResultFh}
#
# RETURN NULL
############################################################################################################################################################     
sub processResults {

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Processing model application results");

    return 0 if ($self->pipelineHasErrors());

    my $model = $self->getModel();

    #Assign scores generated by Model to all peptides. Add both raw scores and benchmark scores (the FPR and TPR for each raw score as it would appear in the training set) 
    my $applicationSetIdList = $self->{ApplicationSetIdList};
    my $resultSet = $model->parseResults($applicationSetIdList);
    return 0 unless $resultSet;
    
    my $foundAllScores = $model->addResultsToPeptideInfo($resultSet, $PROTEOME);
    return 0 unless $foundAllScores;
    
    $model->addBenchmarkScoresToPeptideInfo($PROTEOME);

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "processResults");

}

sub writeResults{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Writing ApplicationPipeline results");

    my $model = $self->getModel();

    #Write results
    my $resultFh = $self->{ResultFh};

    my $columnMap;
    my $lineCounter = 0;

    #write column headers
   
    my $columnHeaderString = $self->makeColumnHeaderString();
    
    print $resultFh $columnHeaderString . "\n";

    my $sequenceCount = 0;
    my $peptideCount = 0;
    
    my $errors = $self->pipelineHasErrors();
    if ($errors){
	$self->writeNoOutput($errors, $columnHeaderString);
    }
    else {
	
	#use ColumnInfo to prepare order that features will be written to the results file
	my $peptides = $model->getPeptides($PROTEOME);
	my $columns = $self->{ColumnInfo};
	
	foreach my $peptideSequenceId (keys %$peptides){
	    $sequenceCount++;
	    my $startPositions = $peptides->{$peptideSequenceId};
	    
	    foreach my $startPosition (keys %$startPositions){
		$peptideCount++;
		my $seqInfo = $startPositions->{$startPosition};
		my $outputLine = $self->makeOutputLine($seqInfo, $columns);
		print $resultFh $outputLine . "\n"; 
	    }
	}
	my $noPeptidesParsedSeqs = $model->getNoPeptidesParsedSeqs();

	#write all seqs with no peptides parsed to output (make temporary hacky $seqInfo hash)
	my $seqsWithNoPeptidesCount = scalar (keys %$noPeptidesParsedSeqs);
	if ($seqsWithNoPeptidesCount > 0){
	    $self->writeLog("$seqsWithNoPeptidesCount sequences had no peptides parsed in the feature pipeline; writing those to output");
	}
	foreach my $modbaseSeqId (keys %$noPeptidesParsedSeqs){
	    my $uniprotAccession = $noPeptidesParsedSeqs->{$modbaseSeqId};
	    my $seqInfo;
	    my $noPeptidesParsedKeyword = $self->getParam("keyword_no_peptides_parsed");

	    $seqInfo->{"Sequence ID"} = $modbaseSeqId;
	    $seqInfo->{"Uniprot Accession"} = $uniprotAccession;
	    $seqInfo->{"Errors"} = $noPeptidesParsedKeyword;

	    my $outputLine = $self->makeOutputLine($seqInfo, $columns);
	    print $resultFh $outputLine . "\n";
	}
    }
    $self->writeLog("Wrote results for $sequenceCount sequences containing $peptideCount total peptides");
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeResults");
}

sub makeOutputLine{
    my ($self, $seqInfo, $columns) = @_; 
    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");
    my $keywordApplicationDisplay = $self->getParam("keyword_application_display");

    my $outputLine = "";
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   
	if ($modes =~ /$keywordApplicationDisplay/){
			
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
    return $outputLine;
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

#############################################################################################################################################################  
# makeColumnHeaderString
# Reads $self->{ColumnInfo} and gets all column names that will be output in the results file. Creates a tab-delimited string with the names. Columns are ordered
# according to specification in $self->getParam("column_info_file") which was loaded previously.
#
# RETURN $columnHeaderString: tab-delimited string with all names of columns to be output.
#############################################################################################################################################################  
sub makeColumnHeaderString{

    my ($self) = @_;

    my $columnHeaderString = "";
    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");
    my $keywordApplicationDisplay = $self->getParam("keyword_application_display");
    
    my $columns = $self->{ColumnInfo};
    
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   
	if ($modes =~ /$keywordApplicationDisplay/){
		    
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

#############################################################################################################################################################  
# writeGlobalError
#
#############################################################################################################################################################  
sub writeGlobalError{
    my ($self, $errorCode, $methodName, $errorMessage) = @_;
    my $errorDescription = $self->{ErrorCodes}->{$errorCode};

    unless($errorDescription){
	my $errorCodes = $self->{ErrorCodes};
	my $errorList = join (", ", keys %$errorCodes);
	my $msg = "Attempted to write an error code ($errorCode) that is not a defined error code. List of allowed codes: $errorList";
	
	$self->writeGlobalError("internal_error", $methodName, $msg . " (original error: $errorMessage)");
	#note this results in an infinite loop if "internal_error" isn't set. There is a check to ensure "internal_error" is set upon initial call to loadErrorCodes()
	
	return 0;  
	#Note that since writeError is called from so many places, there is a good chance that processing will continue in that method. This should be fine although I haven't gone through
	#and tested to make sure that we get expected output if someone writes a bad error code. After that method ends, pipeline will skip to writeUserResults() due to having a global error
	#set here and will proceed with normal global error processing 
    }
  
    my $columnHeaderString = $self->makeColumnHeaderString();
    $self->writeResultFhGlobalError($errorCode, $columnHeaderString);
}





return 1;

