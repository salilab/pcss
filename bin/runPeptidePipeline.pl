#!/usr/bin/perl

use strict;
use FileHandle;

use PeptidePipeline;
use Pipeline;

use Getopt::Long;

my $parameterFileName;

#revert to previous mode -- no iterations


my $options = GetOptions("parameterFileName=s" => \$parameterFileName,
    );

unless ($parameterFileName){
    print STDERR "Usage:  perl runPeptidePipeline.pl --parameterFileName \n";
    exit(1);
}

my $pipeline = PeptidePipeline->new($parameterFileName);

my $serverMode = $pipeline->getParam("server_mode");


if ($serverMode eq "training"){
    $pipeline->readPeptideInputFile();
}
elsif ($serverMode eq "application"){
    
    my $applicationSpec = $pipeline->getParam("application_specification");
    
    if ($applicationSpec eq "S"){
	$pipeline->parsePeptidesFromSequences();

    }
    elsif ($applicationSpec eq "D") {
	$pipeline->readPeptideInputFile();
    }
    #else die!
}

$pipeline->getProteinNames();
    
$pipeline->getBestModels();

$pipeline->parseDsspResults();

$pipeline->runDisopred();

$pipeline->runPsipred();

$pipeline->printAllPeptides();

$pipeline->writeUserResults();


$pipeline->finalize($parameterFileName);

