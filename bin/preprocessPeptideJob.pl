#!/usr/bin/perl

use strict;
use FileHandle;
use FastaReader;

my $jobDirectory = $ARGV[0];
my $parameterFileName = $ARGV[1];

#This should be converted to python to go directly into the job subclass later
unless ($jobDirectory && $parameterFileName){
    print STDERR "preprocessPeptideJob.pl: Exiting; usage: perl preprocessPeptideJob.pl <jobDirectory> <parameterFileName>\n";
    exit(1);
}

my $fullParameterFileName = $jobDirectory . "/" . $parameterFileName;
my $parameters = &readParameterFile($fullParameterFileName);

my $inputPeptideFile = &getParam($parameters, "input_fasta_file_name");
my $sequencesDir = &getParam($parameters, "top_level_seq_directory");
my $fullSequencesDir = $jobDirectory . "/" . $sequencesDir;

my $fullInputPeptideFile = $jobDirectory . "/" . $inputPeptideFile;

my $reader = FastaReader->new($fullInputPeptideFile, 1);
&handleReaderError($fullInputPeptideFile) unless $reader;

$reader->read(); 
my $sequences = $reader->getSequences();

$parameters = &setParam("head_node_preprocess_directory", $jobDirectory, $parameters); #TODO -- change this once we can get the correct job name as a parameter (see Benchmarker.pm)

&writeParameterFile($fullParameterFileName, $parameters);  #overwrite parameter file to add the head_node_preprocess_directory

&writeClusterInput($sequences, $fullSequencesDir, $parameters);

sub writeClusterInput{
    my ($sequenceInfo, $fullSequenceDir, $parameters) = @_;

    my $sequenceCount = 0;

    foreach my $sequenceHeader (keys %$sequenceInfo){
	
	my $residueSequenceList = $sequenceInfo->{$sequenceHeader};
	my $sequence = $residueSequenceList->[0];
	my @cols = split('\|', $sequenceHeader);
	my $modbaseSeqId = $cols[0];

	my $seqNotFoundKeyword = &getParam($parameters, "keyword_no_modbase_for_uniprot");  #TODO -- remove this if we actually aren't passing these sequences from frontend
	next if ($modbaseSeqId eq $seqNotFoundKeyword);

	#write fasta file
	my $outputDir = &makeThreeLetterOutputDir($sequencesDir, $modbaseSeqId);
	my $fastaFileName = &getParam($parameters, "input_fasta_file_name");
	my $fullFastaFileName = $outputDir . "/" .  $fastaFileName;

	my $outputFh = FileHandle->new(">" . $fullFastaFileName);
	unless ($outputFh){
	    print STDERR "preprocessPeptideJob.pl: Exiting; could not open single sequence input file $fullFastaFileName for writing: $!\n";
	    exit(1);
	}
	print $outputFh ">" . $sequenceHeader . "\n" . $sequence . "\n";
	$outputFh->close();

	$sequenceCount++;

	#write parameter file
	my $jobParameterFileName = &getParam($parameters, "job_parameter_file_name");
	my $fullJobParameterFileName = $outputDir . "/" . $jobParameterFileName;
	&writeClusterParameterFile($fullJobParameterFileName, $parameters, $modbaseSeqId);

	print STDERR "preprocessPeptideJob: Wrote three-letter directories in preparation for cluster input; $sequenceCount total sequences are included\n";
    }
} 

sub writeClusterParameterFile{

    my ($jobParameterFileName, $parameters, $modbaseSeqId) = @_;
    $parameters = &setParam("run_name", $modbaseSeqId, $parameters);
    
    &writeParameterFile($jobParameterFileName, $parameters);

}

sub writeParameterFile{   #this is factored out from writeClusterParameterFile in case preprocess needs to rewrite parameter file after adding to it
    my ($parameterFileName, $parameters) = @_;
    my $paramFh = FileHandle->new(">" . $parameterFileName);
    unless ($paramFh){
	print STDERR "preprocessPeptideJob.pl: Exiting; could not open single sequence job parameter file name $parameterFileName for writing: $!\n";
	exit(1);
    }
    
    foreach my $parameterName (keys %$parameters){
	my $parameterValue = $parameters->{$parameterName};
	print $paramFh "$parameterName\t$parameterValue\n";
    }
    $paramFh->close();
}


sub makeThreeLetterOutputDir{
    my ($sequencesDir, $modbaseSeqId) = @_;
    $modbaseSeqId =~ /^(\S\S\S).*/;
    my $threeLetterDir = $1;
    
    my $outputDir = $sequencesDir . "/" . $threeLetterDir . "/" . $modbaseSeqId;
    my $cmd = "mkdir -p $outputDir";
    system($cmd);
    return $outputDir;
}

sub readParameterFile{
    my ($parameterFileName) = @_;

    my $parameterFh = FileHandle->new("<" . $parameterFileName);
    unless ($parameterFh){
	print STDERR "preprocessPeptideJob.pl: Exiting; could not open parameter file name $parameterFileName for reading: $!\n";
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

sub setParam{

    my ($paramName, $paramValue, $parameters) = @_;

    $parameters->{$paramName} = $paramValue;
    #This originally checked to make sure the parameter name didn't already exist, but that caused trouble with adding multiple modbaseSeqIds to the hash
    #since we are setting only three parameters currently, it seems OK to disable the check.
    #otherwise we'd have to create a set of parameters when we first start the script and copy all of
    #them into a secondary parameter hash specific to each modbaseSeqId, which seems like too much
   
    return $parameters;
}


sub handleReaderError{
    my ($inputFile) = @_;
    my $inputFh = FileHandle->new("<" . $inputFile);
    unless ($inputFh){
	print STDERR "preprocessPeptideJob.pl: Exiting due to FastaReader error; could not open fasta file $inputFile for reading: $!\n";
	exit(1);
    }
    else {
	print STDERR "preprocessPeptideJob.pl: Exiting due to FastaReader error; Reader did not return true value when reading $inputFile\n";
	#this should never happen, extreme sanity check
	exit(1);
    }
    
}

sub getParam{

    my ($parameters, $paramName) = @_;
    my $paramValue = $parameters->{$paramName};

    if (!($paramValue)){
        my $errorString = "preprocessPeptideJob.pl: Exiting; tried to retrieve value for parameter $paramName but this is not a valid parameter.  Valid parameters:\n";

        foreach my $parameter (keys %$parameters){
            $errorString .= "--" . $parameter . "\n";
        }
        print STDERR $errorString . "\n";
	exit(1);
    }
    return $paramValue;
}
