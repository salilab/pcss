#!/usr/bin/perl

use strict;
use FileHandle;

my $jobDirectory = $ARGV[0];
my $parameterFileName = $ARGV[1];

unless ($jobDirectory && $parameterFileName){
    print STDERR "postprocessTrainingJob.pl: Exiting; usage: perl postprocessTrainingJob.pl <jobDirectory> <parameterFileName>\n";
    exit(1);
}


my $fullParameterFileName = $jobDirectory . "/" . $parameterFileName;
my $parameters = &readParameterFile($fullParameterFileName);

my $svmResultDir = &getParam($parameters, "top_level_svm_directory");
my $iterationCount = &getParam($parameters, "iteration_count");

my $referencePositiveCount = 0;
my $referenceNegativeCount = 0;
my $iterationNegativeCount;
my @criticalEvalues;
my @criticalTpRates;
my @criticalFpRates;
my $fpsAtEachTp;

&copyUserModelFile($jobDirectory, $parameters, $svmResultDir);

#read all SVM results files to get frequencies of fp rates at all tp rates
for (my $i = 1; $i <= $iterationCount; $i++){
    my $fullSvmResultDir = $jobDirectory . "/" . $svmResultDir . "/svm_iteration_" . $i;  #TODO -- consider parameterizing (think it's OK though)
    my $svmResultFile = &getParam($parameters, "model_pipeline_result_file_name");
    my $fullSvmResultFile = $fullSvmResultDir . "/" . $svmResultFile;
    my $resultFh = FileHandle->new("<" . $fullSvmResultFile);
    unless ($resultFh){
	print STDERR "postprocessPeptideJob.pl: Exiting; could not open svm result file name $fullSvmResultFile for reading\n";
	exit(1);
    }
    my @tpCounts;
    my @fpCounts;

    while (<$resultFh>){
	chomp;
	my $line = $_;
	my @cols = split('\t', $line);
	if (scalar(@cols) == 3){   #tpCount, fpCount, #score
	    my $tpCount = $cols[0];
	    my $fpCount = $cols[1];
	    my $score = $cols[2];
	    $fpsAtEachTp->{$tpCount}->{fpCounts}->{$fpCount}++;
	    $fpsAtEachTp->{$tpCount}->{scores}->{$score}++;
	    
	    push (@tpCounts, $tpCount);
	    push (@fpCounts, $fpCount);

	}
	elsif (scalar(@cols) == 4){  #critical fp rate, critical tp rate, critical evalue, total negatives
	    my $criticalFpRate = $cols[0];
	    my $criticalTpRate = $cols[1];
	    my $criticalEvalue = $cols[2];
	    $iterationNegativeCount = $cols[3];
	    push (@criticalEvalues, $criticalEvalue);
	    push (@criticalTpRates, $criticalTpRate);
	    push (@criticalFpRates, $criticalFpRate);
	    
	}
	else {
	    my $errorColumnCount = scalar(@cols);
	    print STDERR "postprocessPeptideJob.pl: Exiting; got unexpected number of columns ($errorColumnCount) when reading svm training result file $fullSvmResultFile (expect 2 or 3 columns)\n";
	    exit(1);
	}
    }
    my $countSize = scalar(@tpCounts);
    my $positiveCount = $tpCounts[$countSize - 1];   #max value of @tpCounts is the number of positives in the test set
    my $negativeCount = $fpCounts[$countSize - 1];

    $referencePositiveCount = &checkReferenceTestSetCount($referencePositiveCount, $positiveCount, "positive"); 
    $referenceNegativeCount = &checkReferenceTestSetCount($referenceNegativeCount, $iterationNegativeCount, "negative"); 
}

#now go through all FP Counts at each TP count to get average rates
my $finalResultFileName = &getParam($parameters, "training_final_result_file_name");  #TODO -- see if it's a better idea to have all output file be the same name regardless of server mode
my $fullFinalResultFileName = $jobDirectory . "/" . $finalResultFileName;
my $finalResultFh = FileHandle->new(">" . $fullFinalResultFileName);
unless ($finalResultFh){
    print STDERR "postprocessPeptideJob.pl: Exiting; could not open final training result file name $fullFinalResultFileName for writing results\n";
    exit(1);
}

print $finalResultFh "fpr\ttpr\n";
print $finalResultFh "0.0\t0.0\n";

my $totalFps = $referenceNegativeCount * $iterationCount;

foreach my $tp (sort ({$a <=> $b} keys %$fpsAtEachTp)){
    my $fpCounts = $fpsAtEachTp->{$tp}->{fpCounts};
    my $scores = $fpsAtEachTp->{$tp}->{scores};
    my $totalFpsAtThisTp = 0;
    
    my @allFps;
	
    foreach my $fpCount (keys %$fpCounts){
	my $fpCountValue = $fpCounts->{$fpCount};
	$fpCount *= $fpCountValue;
	$totalFpsAtThisTp += $fpCount;

	my $fpRate = ($fpCount * 1.0)  / ($totalFps * 1.0);  
	push (@allFps, $fpRate);  #prepare standard deviation list
    }
    my $fpRateAtThisTp = ($totalFpsAtThisTp * 1.0) / ($totalFps * 1.0);
    my $tpRate = ($tp * 1.0) / ($referencePositiveCount * 1.0);


    #get standard deviation
    my $fpRateSum = 0;
    foreach my $fpRate (@allFps){
	$fpRateSum += $fpRate;
    }
    my $count = scalar(@allFps);
    my $fpRateAverage = ($fpRateSum * 1.0) / ($count * 1.0);

    my $averageDifferenceSum = 0;
    print STDERR "fpRateAverage: $fpRateAverage\n";
    foreach my $fpRate (@allFps){
	print STDERR "next fp rate: $fpRate\n";
	my $difference = $fpRate - $fpRateAverage;
	my $squaredDifference = $difference * $difference;
	$averageDifferenceSum += $squaredDifference;
    }
    print STDERR "final average difference sum: $averageDifferenceSum count $count\n";
    $averageDifferenceSum /= $count;
    my $stdDev = sqrt($averageDifferenceSum);

    #get average score at this tp
    my $scoreTotal = 0.0;
    
    foreach my $score (keys %$scores){
	$scoreTotal += $score;
    }
    my $scoreAverage = $scoreTotal / ($iterationCount * 1.0);
    

    print $finalResultFh "$fpRateAtThisTp\t$tpRate\t$scoreAverage\t$stdDev\n";    

}
       
print $finalResultFh "1.0\t1.0\n\n\n";

#print $finalResultFh "***False Positive Rates, True Positive Rates at each SVM score:\n";

my $totalEvalue = 0;
foreach my $eValue (@criticalEvalues){
    
    $totalEvalue += $eValue;
}
my $avgEvalue = $totalEvalue / ($iterationCount * 1.0);
print $finalResultFh "\naverage critical evalue: $avgEvalue\n";


my $totalTpAtEvalue = 0;
foreach my $tpRate (@criticalTpRates){
    $totalTpAtEvalue += $tpRate;
}

my $avgTpAtEvalue = $totalTpAtEvalue / ($iterationCount * 1.0);
print $finalResultFh "\naverage tp rate at critical eValue: $avgTpAtEvalue\n";

my $totalFpAtEvalue = 0;
foreach my $fpRate (@criticalFpRates){
    $totalFpAtEvalue += $fpRate;
}
my $avgFpAtEvalue = $totalFpAtEvalue / ($iterationCount * 1.0);
print $finalResultFh "\naverage fp rate at critical eValue: $avgFpAtEvalue\n";



sub checkReferenceTestSetCount{

    #sanity check to make sure each SVM iteration has the same number of peptides that it tested on
    my ($currentReferenceCount, $countToCheck, $tag) = @_;
    if ($currentReferenceCount == 0){
	$currentReferenceCount = $countToCheck;
    }
    else {
	if ($countToCheck != $currentReferenceCount){
	    print STDERR "postprocessPeptideJob.pl: Exiting; did not have consistent number of $tag peptides across different svm iterations (found one with $countToCheck and another with $currentReferenceCount\n";
	    exit(1);
	}
    }
    return $currentReferenceCount;

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

sub copyUserModelFile{
    my ($jobDirectory, $parameters, $svmResultDir) = @_;

    my $svmModelFileName = &getParam($parameters, "svm_complete_model_name");
    my $fullSvmModelFile = $jobDirectory . "/" . $svmResultDir . "/svm_iteration_1/" . $svmModelFileName;  #TODO - parameterize this; sync with __init__.py
    my $cpCmd = "cp $fullSvmModelFile $jobDirectory";
    print STDERR "copying user model file with $cpCmd\n";
    system($cpCmd);
}


#   foreach my $fpRate (@allFps){
#	
#	my $squaredFp = $fpRate * $fpRate;
#	$squaredFpSum += $squaredFp;
#   }
#   my $count = scalar(@allFps);
#   my $tempFpAverage = ($squaredFpSum * 1.0) / ($count * 1.0);
#   $squaredFpSum /= $count;
#   my $squaredAverage = $tempFpAverage * $tempFpAverage;
#   my $stdDev = sqrt($squaredFpSum - $squaredAverage);#

    #get average score at this tp
#   my $scoreTotal = 0.0;
    
#   foreach my $score (keys %$scores){
#	$scoreTotal += $score;
#   }
#   my $scoreAverage = $scoreTotal / ($iterationCount * 1.0);
    

   #print $finalResultFh "$fpRateAtThisTp\t$tpRate\t$scoreAverage\t$stdDev\n";
