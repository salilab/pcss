package Pipeline;

use BenchmarkerPipeline;
#use CreationPipeline;
use ApplicationPipeline;
use PeptidePipeline;
use POSIX qw(ceil floor);

use strict;
use FileHandle;

sub new{
    my ($class, $parameterFileName, $type) = @_;
    my $self = {};
    bless $self, $class;

    #read params
    $self->readParameterFile($parameterFileName);

    #set global pipeline data structures
    $self->loadLog($type);
    $self->loadResultFh($type);
    $self->loadColumnInfo();

    $self->{PipelineType} = $type; #todo- check

    return $self;
}


sub readParameterFile{

    my ($self, $parameterFile) = @_;
    my $parameterFh = FileHandle->new("<" . $parameterFile) || die "could not open parameterFile $parameterFile\n";

    my $parameters;

    while (<$parameterFh>){
	chomp;
	my $line = $_;
	next if ($line =~ /^\#/);

	my ($paramName, $paramValue) = split('\t', $line);
	$parameters->{$paramName} = $paramValue;
    }
    $self->setParams($parameters);
}


sub setParams{
    my ($self, $parameters) = @_;
    $self->{Parameters} = $parameters;
}


sub getParams{
    my ($self) = @_;
    my $params = $self->{Parameters};
    return $params;
}


sub loadPeptides{
    my ($self) = @_;
    print STDERR "method loadPeptides() should only be called from subclass of Pipeline\n";
    exit(1);
}

sub execute{
    my ($self) = @_;
    print STDERR "method execute() should only be called from subclass of Pipeline\n";
    exit(1);
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


sub finalize{

    my ($self, $parameterFileName) = @_;
    
    my $logFh = $self->{Log};
    $logFh->close();
    my $resultFh = $self->{ResultFh};
    $resultFh->close();
 
}


sub getParam{

    my ($self, $paramName) = @_;
    my $paramValue = $self->{Parameters}->{$paramName};
    if (!($paramValue)){
 	print STDERR "Error: given parameter ($paramName) is not valid.  Valid parameters";

	my $logString = "";
	foreach my $parameter (keys %{$self->{Parameters}}){
	    $logString .= "--" . $parameter . "\n";
	}
	print STDERR "$logString\n";

	exit(1);
    }
    return $paramValue;
}


sub loadLog{
    my ($self, $type, $runName) = @_;

    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $logFileParamName = $type . "_pipeline_log_name";
    my $logFileName =$self->getParam($logFileParamName);
    my $fullLogFileName = $pipelineDir . "/" . $run . "/" . $logFileName;

    my $logFh = FileHandle->new(">" . $fullLogFileName) || die "could not open $fullLogFileName for writing\n";
    $self->{LogFileName} = $fullLogFileName;
    $self->{Log} = $logFh;
}


sub loadResultFh{
    my ($self, $type, $runName) = @_;

    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $resultsFileParamName = $type . "_pipeline_result_file_name";
    my $resultsFileName = $self->getParam($resultsFileParamName);
    my $fullResultsFileName = $pipelineDir . "/" . $run . "/" . $resultsFileName;
    my $resultsFh = FileHandle->new(">>" . $fullResultsFileName) || die "could not open $fullResultsFileName for writing\n";
    $self->{ResultsFileName} = $fullResultsFileName;
    $self->{ResultFh} = $resultsFh;
}


sub writeResultFhOutputError{
    my ($self) = @_;
    my $resultFh = $self->{ResultFh};
    $resultFh->close();

    my $type = $self->{PipelineType};
    
    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $resultsFileParamName = $type . "_pipeline_result_file_name";   #todo -- test this for different pipeline types; not sure if it's entirely accuratey
    my $resultsFileName = $self->getParam($resultsFileParamName);
    my $fullResultsFileName = $pipelineDir . "/" . $run . "/" . $resultsFileName;
    my $resultsFh = FileHandle->new(">" . $fullResultsFileName) || die "could not open $fullResultsFileName for writing\n";

    my $outputErrorKeyword = $self->getParam("keyword_output_error");

    print $resultsFh $outputErrorKeyword;
}

sub writeResultFhGlobalError{
    my ($self, $errorName, $columnHeaderString) = @_;
    my $resultFh = $self->{ResultFh};
    $resultFh->close();

    my $resultFileName = $self->{ResultsFileName};
    $resultFh = FileHandle->new(">" . $resultFileName) || die "could not open $resultFileName for writing\n";

    print $resultFh $columnHeaderString . "\n" . $errorName;

}
    


sub writeLog{
    my ($self, $message) = @_;
    my $logFh = $self->{Log};
    my ($sec,$min,$hour,$mday,$mon,$year,$wday, $yday,$isdst)=localtime(time);

    my $dateLine = sprintf ("%4d-%02d-%02d %02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

    my $cutPoint = 90;
    my $previousCutPoint = 0;
    my @lines;
    while (1){
        my $currentCharAt = substr($message, $cutPoint, 1);
        while ($currentCharAt =~ /\S/){
            $cutPoint++;
            $currentCharAt = substr($message, $cutPoint, 1);
        }

        my $substr = substr($message, $previousCutPoint, $cutPoint - $previousCutPoint);
        $previousCutPoint = $cutPoint;
        $cutPoint += 90;
        push(@lines, $substr);
        last if ($previousCutPoint > length($message));
    }
    my $firstLine = $lines[0];
    print $logFh "$dateLine:\t $firstLine\n";
    for(my $i = 1; $i < scalar(@lines); $i++){
        my $line = $lines[$i];
        print $logFh  "\t\t\t$line\n";
    }
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
