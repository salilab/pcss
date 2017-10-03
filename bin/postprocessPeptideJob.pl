#!/usr/bin/perl

use strict;
use FileHandle;
use FastaReader;


#1. Go to base data dir where all files have been copied.
#2. Copy each line by line into global output file if appropriate according to column info parameter file

my $jobDirectory = $ARGV[0];
my $parameterFileName = $ARGV[1];

#This should be converted to python to go directly into the job subclass
unless ($jobDirectory && $parameterFileName){
    print STDERR "postprocessPeptideJob.pl: Exiting; usage: perl preprocessPeptideJob.pl <jobDirectory> <parameterFileName>\n";
    exit(1);
}

#read parameters
my $fullParameterFileName = $jobDirectory . "/" . $parameterFileName;
my $parameters = &readParameterFile($fullParameterFileName);

#get input fasta file and other params
my $seqDirectory = &getParam($parameters, "top_level_seq_directory");
my $baseDataDir = $jobDirectory . "/" . $seqDirectory;
my $inputFastaFile = &getParam($parameters, "input_fasta_file_name");
my $fullInputFastaFile = $jobDirectory . "/" . $inputFastaFile;
my $peptidePipelineResultsFile = &getParam($parameters, "model_pipeline_result_file_name");
my $noPeptidesParsedKeyword = &getParam($parameters, "keyword_no_peptides_parsed");

#open final output file
my $applicationFinalResultFileName = &getParam($parameters, "application_final_result_file_name");
my $fullAppFinalResultFileName = $jobDirectory . "/" . $applicationFinalResultFileName;
my $finalOutputFh = FileHandle->new(">" . $fullAppFinalResultFileName);
unless ($finalOutputFh){
    print STDERR "postprocessPeptideJob.pl: Exiting; could not open final application output file $fullAppFinalResultFileName for writing\n";
    exit(1);
}


#open temp logging file
my $tempLogFile = $jobDirectory . "/postProcessLog.txt";
my $tempLogFh = FileHandle->new(">" . $tempLogFile);
unless ($tempLogFh){
    print STDERR "postprocessPeptideJob.pl: Exiting; could not open temp log file $tempLogFile for writing\n";
    exit(1);
}

#read in input file to get sequence IDs
require FastaReader;
print STDERR "reading sequences from $inputFastaFile\n";
my $reader = FastaReader->new($fullInputFastaFile, 1);
&handleReaderError($fullInputFastaFile) unless $reader;
$reader->read();
my $sequences = $reader->getSequences();

#read in columns from param file; write display names to output
my $columnInfoFile = &getParam($parameters, "column_info_file");
my $keywordApplicationDisplay  = &getParam($parameters, "keyword_application_display");  # TODO - change this if there are other column modes besides those server modes
my $columnInfo = &getColumnInfo($columnInfoFile);
my $columnHeaderString = &makeColumnHeaderString($columnInfo, $keywordApplicationDisplay, $parameters);  
print $finalOutputFh $columnHeaderString . "\n";

my $seqCount = 0;

#go through all input sequences to get application pipeline results for each
foreach my $sequenceHeader (keys %$sequences){
    $seqCount++;
    if ($seqCount % 1000 == 0){
	my $time = localtime();
	print $tempLogFh "$time: processing sequence $seqCount\n";
    }
    #use sequence name and first three letters to create output directory
    my ($modbaseSeqId, $uniprotAccession) = split('\|', $sequenceHeader);
    $modbaseSeqId =~ /^(\S\S\S).*/;
    my $threeLetterDir = $1;
    my $outputDir = $baseDataDir . "/" . $threeLetterDir . "/" . $modbaseSeqId;

    #get filehandle to results for this sequence
    my $fullResultsFile = $outputDir . "/" . $peptidePipelineResultsFile;    
    unless (-e $fullResultsFile){
	print $finalOutputFh $uniprotAccession . "\tThis accession could not be processed due to an internal error.  A report has been generated and the service developers have been notified.  We apologize for the inconvenience (modbase seq id $modbaseSeqId)\n";
	next;
    }
    my $resultsFh = FileHandle->new("<" . $fullResultsFile);
    unless ($resultsFh){
	print STDERR "postprocessPeptideJob.pl: Exiting; could not open individual sequence results file $fullResultsFile\n";
	exit(1);
    }
    
    #read first line of result file to get column headers and map them to their order in the result file.  
    #TODO -- I think this can be avoided by looking at global column info and assuming that it will match up with each result file, but need to be sure
    my $sequenceFileColumnLine = <$resultsFh>;
    my $columnNumberToNameMap;
    my $counter = 1;
    my @sequenceFileColumnNames = split('\t', $sequenceFileColumnLine);
    foreach my $sequenceFileColumnName (@sequenceFileColumnNames){
	$columnNumberToNameMap->{$counter} = $sequenceFileColumnName;
	$counter++;
    }
     
    #go through actual results
    while (<$resultsFh>){
	chomp;
	my $line = $_;

	#check there were peptides to process
	if ($line =~ /$noPeptidesParsedKeyword/){
	    print $finalOutputFh "No peptides were parsed from uniprot accession $uniprotAccession given the rules read from the user-provided peptide specifier file\n";
	    last;
	}
	
	#read result line for peptide, write to $outputInfo (hash of the form $outputInfo->{$columnName} = $value
	my @values = split('\t', $line);
        $counter = 1;
	my $outputInfo;
	foreach my $value (@values){  
	    my $columnNameForThisValue = $columnNumberToNameMap->{$counter};   #get column name for the current value
	    $outputInfo->{$columnNameForThisValue} = $value;  
	    $counter++;
	}

	#create and write line for this peptide
	my $finalLine = &makeOutputline($outputInfo, $columnInfo, $keywordApplicationDisplay, $parameters);
	print $finalOutputFh $finalLine . "\n";
    }
}


sub makeOutputline{
    
    my ($outputInfo, $columns, $peptideColumnModeName, $parameters) = @_;

    my $outputLine = "";
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   #check if the output file for Peptide Pipeline mode is supposed to output this column type

	if ($modes =~ /$peptideColumnModeName/){
	    my $method = $columns->{$columnShortName}->{method};
	    my $columnName = $columns->{$columnShortName}->{displayName};
	    if ($method eq "universal"){
	    
		my $nextOutputValue = $outputInfo->{$columnName};
		$outputLine .= $nextOutputValue . "\t";
	    }
	    else {
		if (&runMethod($method, $parameters)){
		    my $nextOutputValue = $outputInfo->{$columnName};
		    $outputLine .= $nextOutputValue . "\t";
		}
	    }
	}
    }
    return $outputLine;
}

sub getColumnInfo{

    my ($columnInfoFile) = @_;
    my $fh = FileHandle->new("<" . $columnInfoFile);
    unless ($fh){
	print STDERR "postprocessPeptideJob.pl: Exiting; could not open individual column info file $columnInfoFile\n";
	exit(1);
    }
    
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
    return $columnInfo;

}

sub makeColumnHeaderString{

    my ($columns, $peptideColumnModeName, $parameters) = @_;
   
    my $columnHeaderString = "";

    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   #check if the output file for Peptide Pipeline mode is supposed to output this column type
	if ($modes =~ /$peptideColumnModeName/){
	    my $displayName = $columns->{$columnShortName}->{displayName};
	    my $method = $columns->{$columnShortName}->{method};
	    if ($method eq "universal"){
		$columnHeaderString .= $displayName . "\t";
	    }
	    else {
		if (&runMethod($method, $parameters)){
		    $columnHeaderString .= $displayName . "\t";
		}
	    }
	}
    }
    return $columnHeaderString;
}


sub runMethod{

    my ($methodName, $parameters) = @_;
    
    my $runMethodParamValue = &getParam($parameters, $methodName);

    if ($runMethodParamValue eq "yes"){
	return 1;
    }
    else {
	return 0;
    }
}



sub readParameterFile{
    my ($parameterFileName) = @_;

    my $parameterFh = FileHandle->new("<" . $parameterFileName);
    unless ($parameterFh){
	print STDERR "postprocessPeptideJob.pl: Exiting; could not open parameter file name $parameterFileName for reading\n";
	exit(1);
    }
    
    my $parameters;

    while (<$parameterFh>){
        chomp;
        my $line = $_;
        my ($paramName, $paramValue) = split('\t', $line);
	$parameters->{$paramName} = $paramValue;
    }

    return $parameters;
}

sub getParam{

    my ($parameters, $paramName) = @_;
    my $paramValue = $parameters->{$paramName};

    if (!($paramValue)){
        my $errorString = "postprocessPeptideJob.pl: Exiting; tried to retrieve value for parameter $paramName but this is not a valid parameter.  Valid parameters:\n";

        foreach my $parameter (keys %$parameters){
            $errorString .= "--" . $parameter . "\n";
        }
        print STDERR $errorString . "\n";
	exit(1);
    }
    return $paramValue;
}

sub handleReaderError{
    my ($inputFile) = @_;
    my $inputFh = FileHandle->new("<" . $inputFile);
    unless ($inputFh){
	print STDERR "postprocessPeptideJob.pl: Exiting due to FastaReader error; could not open fasta file $inputFile for reading: $!\n";
	exit(1);
    }
    else {
	print STDERR "postprocessPeptideJob.pl: Exiting due to FastaReader error; Reader did not return true value when reading $inputFile\n";
	#this should never happen, extreme sanity check
	exit(1);
    }
    
}

