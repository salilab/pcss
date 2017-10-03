#!/usr/bin/perl

use strict;
use FileHandle;


#Takes results of "select * from newmodels where run='human_2008'\G which I saved to a file; parses them and makes a tab-delimited table.
my $fileName = $ARGV[0];

unless ($fileName){
    print STDERR "usage: perl convertHuman2008results.pl <fileName>\n";
    exit(1);
}


my $fh = FileHandle->new("<" . $fileName) || die "could not open $fileName\n";

my $resultsHash;
my $firstRow = <$fh>;

my $currentRowNumber;

if ($firstRow =~ /\*\*\*\s(\d+)\.\srow/){
    $currentRowNumber = $1;
}

my $counter = 0;
my $valueLine = "";
while (<$fh>){
    chomp;
    my $line = $_;

    if ($line =~ /\*\*\*\s(\d+)\.\srow/){
	$currentRowNumber = $1;               #key all data on row number   *************************** 1. row ***************************\y
	$counter++;
	print STDERR $valueLine . "\n";
	$valueLine = "";
    }
    
    else {
	if ($line =~ /\s*(.*):\s(.*)\\/){      #  align_id: d22a07fbb673b35abc675792572642ec\
	    my $value = $2;
	    $value =~ s/^\s+//; 
	    $value =~ s/\s+$//;
	    $valueLine .= "$value\t";
	}
    }
}
