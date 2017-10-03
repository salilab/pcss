package CleavageSiteModel;

#use HmmModel;
use SvmModel;
use strict;
use FileHandle;

sub new{
    my ($class) = @_;
    my $self = {};
    bless $self, $class;

    return $self;
}


#things to keep in mind.
# -- data structure that holds all of the cleavage site information and its properties
# -- how to handle positives, negatives in training set when not all of them have it.
# -- how to generalize results gathering

sub setParams{

    my ($self, $parameters) = @_;
    $self->{Parameters} = $parameters;
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


sub trainsOnNegatives{
    my ($self) = @_;
    print STDERR "method trainsOnNegatives should only be called from a subclass of CleavageSiteModel\n";

}


#should be "add to peptide count"
sub setPeptideCount{
    my ($self, $count, $type) = @_;

    $self->{PeptideCounts}->{$type} += $count;
}

sub getPeptideCount{

    my ($self, $type) = @_;
    my $count = $self->{PeptideCounts}->{$type};
    return $count;
}

sub trainModel{
    my ($self, $trainingSet, $trainingSetFh) = @_;
    print STDERR "method train() should only be called on a subclass of CleavageSiteModel\n";
    exit(1);
}

sub testModel{

    my ($self, $testSet) = @_; 
    print STDERR "method testModel() should only be called on a subclass of CleavageSiteModel\n";
    exit(1);
}


sub getPeptideSource{
    my ($self, $inputProteinId, $pOnePosition) = @_;
    print STDERR "method getSource() should only be called from subclass of CleavageSiteModel\n";

    exit(1);

}

sub loadPeptides{

    my ($self) = @_;
    print STDERR "method loadPeptides() should only be called from subclass of CleavageSiteModel\n";
    exit(1);
}


sub printPeptides{

    my ($self) = @_;
    print STDERR "method printPeptides() should only be called from subclass of CleavageSiteModel\n";
    exit(1);
}

sub setLog{
    my ($self, $log) = @_;
    $self->{Log} = $log;
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


sub setResultFh{
    my ($self, $resultFh) = @_;
    $self->{ResultFh} = $resultFh;
}



return 1;
