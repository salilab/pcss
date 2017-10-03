package FastaReader;

use strict;
use FileHandle;

#General Usage:
#1. Create new $reader with input file and joinSequences set.
#2. $reader->read();
#3. my $sequences = $reader->getSequences();
#   $sequences is hash with keys as fasta headers and values are arrays containing lines or one element containing whole sequences (see joinSequence below).  No ">" are included

#joinSequence:  if true, value of getSequences() is hash with header as key and one item array containing full sequence.
#               if false, value is array with each entry being one line of the sequence
sub new{

    my ($class, $inputFile, $joinSequence) = @_;

    my $self = {};
    bless $self, $class;

    $self->{inputFile} = $inputFile;
    $self->{joinSequence} = $joinSequence;
    $self->{sequences} = {};

    my $inputFh = FileHandle->new("<" . $inputFile) || return $!;
    $self->{inputFh} = $inputFh;
    return $self;

}

sub read{

    my ($self) = @_;


    my $currentHeaderLine;
    my $sequences;

    my $inputFh = $self->{inputFh};
    while (<$inputFh>){
	chomp;
	$currentHeaderLine = $_;
	if ($currentHeaderLine =~ /\>(.*)/){   #scroll down to first header line; usually is the first line in the file
	    $currentHeaderLine = $1;
	    last;
	}
    }

    my @seqArray;

    #precondition: $currentHeaderLine equal to a header line 
    #(either comes in from previous while loop or end of within this loop with $currentLine set)
    while (<$inputFh>){

	chomp;
	my $line = $_;

	if ($line =~ /^\s*$/){  #skip blank lines
	    next;
	}
	
	if ($line =~ /\>(.*)/){   #reached next entry, end current sequence
	
	    $sequences = $self->addSequence($currentHeaderLine, @seqArray);  #add to list
   	    
	    #reinitialize values
	    $currentHeaderLine = $1;
	    my @newSeqArray;
	    @seqArray = @newSeqArray;

	    next;
	}

	push (@seqArray, $line);
    }

    #process final sequence
    $self->addSequence($currentHeaderLine, @seqArray);
    $inputFh->close();
    return 1;
}


sub getSequences{
    my ($self) = @_;

    my $sequences = $self->{sequences};

    return $sequences;


}


sub addSequence{
    my ($self, $currentHeaderLine, @seqArray) = @_;

    my @finalSeqArray;
    
    if ($self->{joinSequence}){
	my $finalSeqString = join("", @seqArray);
	push (@finalSeqArray, $finalSeqString);
    }
    else {
	@finalSeqArray = @seqArray;
    }
    
    $self->{sequences}->{$currentHeaderLine} = \@finalSeqArray;


}


return 1;
