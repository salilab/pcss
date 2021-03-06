#!/usr/bin/perl

use Pipeline;
use CleavageSiteModel;
use BenchmarkerPipeline;
use ApplicationPipeline;
use CreationPipeline;

use strict;
use Getopt::Long;

my $parameterFileName;
my $pipelineClass;

my $options = GetOptions("parameterFileName=s" => \$parameterFileName,
			 "pipelineClass=s" => \$pipelineClass,
			 
    );


print STDERR "pipeline class: $pipelineClass\n";

unless ($parameterFileName && $pipelineClass){
    print STDERR "Usage:  perl runModelPipeline.pl --parameterFileName --pipelineClass\n";
    print STDERR "Current pipeline classes:  BenchmarkerPipeline, CreationPipeline, ApplicationPipeline\n";
    exit(1);
}

my $pipeline = $pipelineClass->new($parameterFileName);

#$pipeline->loadPeptides();

#$pipeline->printPeptides();

$pipeline->execute();

#$pipeline->readResults();

$pipeline->finalize($parameterFileName);


