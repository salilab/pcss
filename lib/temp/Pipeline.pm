package Pipeline;

use BenchmarkerPipeline;
#use CreationPipeline;
use ApplicationPipeline;
use POSIX qw(ceil floor);

use strict;
use FileHandle;

sub new{
    my ($class, $parameterFileName) = @_;
    my $self = {};
    bless $self, $class;

    $self->initialize($parameterFileName);
    return $self;
}


sub initialize{
    my ($self, $parameterFileName) = @_;

    #read params
    $self->readParameterFile($parameterFileName);

    #set global pipeline data structures
    $self->loadLog();
    $self->loadResultFh();
    my $modelClass = $self->getParam("model_class");

    #create and initialize model object
    my $model = $modelClass->new();
    $self->setModel($model);
    $model->setParams($self->getParams());
    $model->setLog($self->{Log});
    $model->setResultFh($self->{ResultFh});
    $model->initialize();  

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
    my ($self) = @_;

    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $logFileName =$self->getParam("model_pipeline_log_name");
    my $fullLogFileName = $pipelineDir . "/" . $run . "/" . $logFileName;
 
    my $logFh = FileHandle->new(">" . $fullLogFileName) || die "could not open $fullLogFileName for writing\n";
    $self->{LogFileName} = $fullLogFileName;
    $self->{Log} = $logFh;
}


sub loadResultFh{
    my ($self) = @_;

    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $resultsFileName = $self->getParam("model_pipeline_result_file_name");
    my $fullResultsFileName = $pipelineDir . "/" . $run . "/" . $resultsFileName;
    my $resultsFh = FileHandle->new(">>" . $fullResultsFileName) || die "could not open $fullResultsFileName for writing\n";
    $self->{ResultsFileName} = $fullResultsFileName;
    $self->{ResultFh} = $resultsFh;

}

sub writeLog {

    my ($self, $message) = @_;
    my $logFh = $self->{Log};
    my ($sec,$min,$hour,$mday,$mon,$year,$wday, $yday,$isdst)=localtime(time);

    my $dateLine = sprintf ("%4d-%02d-%02d %02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

    print $logFh "$dateLine (Pipeline):\t$message\n";
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
