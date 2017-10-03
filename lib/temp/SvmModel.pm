package SvmModel;

use strict;
use FileHandle;

use CleavageSiteModel;
use POSIX qw(ceil floor);
use List::Util qw(sum);

our @ISA = qw(CleavageSiteModel);

my $POSITIVE = "positive";
my $NEGATIVE = "negative";
my $TEST = "Test";
my $PROTEOME = "Application";


sub new{

    my ($class) = @_;
    my $self = $class->SUPER::new();
    bless ($self, $class);
    
    return $self;
}

sub loadPeptides{

    my ($self, $peptideFileName) = @_;
    
    my $columnMap;    
    my $peptideFh = FileHandle->new("<" . $peptideFileName) || die "could not open $peptideFileName\n";
    my $firstLine = <$peptideFh>;
   
    #Read in column names from first line
    my @cols = split('\t', $firstLine);
    my $i = 0;
    foreach my $col (@cols){
	$self->{ResultFieldStrings}->{$col} = 1;
	$columnMap->{$i} = $col;
	$i++;
    }

    my $noPeptidesParsedKeyword = $self->getParam("keyword_no_peptides_parsed");

    my $resultCounts;
    #read rest of result file
    while (<$peptideFh>){
	
	chomp;
	my $line = $_;
	if ($line =~ /$noPeptidesParsedKeyword/){
	    $self->writeLog("found line $line; writing no peptides parsed");
	    $self->setNoPeptidesParsed();
	    last;
	}
    
	my @cols = split('\t', $line);
	
	my $currentSeqInfo;
	my $i = 0;
	foreach my $colValue (@cols){     #read in field values 
	    my $currentColumn = $columnMap->{$i};
	    $currentSeqInfo->{$currentColumn} = $colValue;
	    $i++;
	}
	my $seqId = $currentSeqInfo->{"Sequence ID"};
	my $startPosition = $currentSeqInfo->{"Peptide Start"};
	my $peptideType = $currentSeqInfo->{"Classification"};
	$resultCounts->{$peptideType}++;
#	$currentSeqInfo->{status} = $peptideType;  #TODO -- see if this is ever used
	$self->{Peptides}->{$peptideType}->{$seqId}->{$startPosition} = $currentSeqInfo;  #will need to test this out
    }
    foreach my $peptideType (keys %$resultCounts){
	my $resultCount = $resultCounts->{$peptideType};
	$self->setPeptideCount($resultCount, $peptideType);
    }

}



sub createPeptideFeatures{
    my ($self, $peptideType) = @_;

    $self->resetFeatureNumber();

    $self->createPeptideSequenceFeatures($peptideType);
    
    $self->createDisorderFeatures($peptideType);

    $self->createStructureAccessibilityFeatures($peptideType);
    
    $self->createSequenceSecondaryStructureFeatures($peptideType);
    
    $self->createStructureSecondaryStructureFeatures($peptideType);
}


sub createPeptideSequenceFeatures{
    
    my ($self, $peptideType) = @_;

    my $peptides = $self->getPeptides($peptideType);
    
    my $maxFeatureNumber = $self->{MaxFeatureNumber};

    my $zeroes = $self->getParam("zeroes");

    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};
	
	foreach my $pOnePosition (keys %$pOnePositions){

	    my $featureString = "";
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $peptideSequence = $self->getResultField($seqInfo, "Peptide Sequence");
	    
	    my @seqArray = split('', $peptideSequence);
    
	    my $counter = 1;
	    my $featureNumber;
	    #Each residue gets feature number corresponding to its identity and position in the peptide sequence.
	    #First position is 1 through 20, second is 21 through 40, etc.
	    foreach my $nextResidue (@seqArray){
		my @residueArray = @{$self->{ResidueStrings}->{$nextResidue}};
				
		foreach my $nextFeature(@residueArray){
		    $featureNumber = $counter + $maxFeatureNumber;
		    $counter++;
		    next if ($zeroes eq "implicit" && $nextFeature == 0);
		    $featureString .= $featureNumber . ":" . $nextFeature . " ";
		}
	    }
	    $seqInfo->{allFeatures}->{peptideSeqFeature} = $featureString;
	}
    }
    my $peptideLength = $self->{PeptideLength};
    $maxFeatureNumber += $peptideLength * 20;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}

sub createDisorderFeatures{

    my ($self, $peptideType) = @_;
    return unless $self->useStructureFeatures();

    my $peptides = $self->getPeptides($peptideType);
    my $maxFeatureNumber = $self->{MaxFeatureNumber};
    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    my $counter = 1;
	    my $featureString = "";
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $disorderString = $self->getResultField($seqInfo, "Disopred Scores");
	    my @disorderScores = split(',', $disorderString);
	    foreach my $disorderScore (@disorderScores){

		my $currentFeatureNumber = $maxFeatureNumber + $counter;
		$featureString .= $currentFeatureNumber . ":$disorderScore ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{disorderFeature} = $featureString;
	}
    }
    my $peptideLength = $self->{PeptideLength};
    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}


sub createSequenceSecondaryStructureFeatures{

    my ($self, $peptideType) = @_;

    return unless $self->useStructureFeatures();

    my $peptides = $self->getPeptides($peptideType);
    
    my $maxFeatureNumber = $self->{MaxFeatureNumber};

    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    my $counter = 1;
	    my $featureString = "";
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $secondaryStructureString = $self->getResultField($seqInfo, "PSIPRED Scores");
	    my @psipredScores = split(',', $secondaryStructureString);
	    foreach my $psipredScore (@psipredScores){

		my $currentFeatureNumber = $maxFeatureNumber + $counter;  
		$featureString .= $currentFeatureNumber . ":$psipredScore ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{sequenceSecondaryStructureFeature} = $featureString;
	}
    }
    my $peptideLength = $self->{PeptideLength};
    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}



sub createStructureSecondaryStructureFeatures{

    my ($self, $peptideType) = @_;
    return unless $self->useStructureFeatures();

    my $peptides = $self->getPeptides($peptideType);
    
    my $maxFeatureNumber = $self->{MaxFeatureNumber};

    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};
 
	foreach my $pOnePosition (keys %$pOnePositions){
	    
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $featureString = "";
	    my $structureSecondaryStructureString = $self->getResultField($seqInfo, "Peptide Structure Values");
	    my $counter = 1;

	    my @secondaryStructureValues = split(',', $structureSecondaryStructureString);
	    foreach my $secondaryStructureValue (@secondaryStructureValues){

		my $currentFeatureNumber = $maxFeatureNumber + $counter;
		$featureString .= $currentFeatureNumber . ":$secondaryStructureValue ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{structureSecondaryStructureFeature} = $featureString;
	}
    }
    my $peptideLength = $self->{PeptideLength};
    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}


sub createStructureAccessibilityFeatures{

    my ($self, $peptideType) = @_;
    return unless $self->useStructureFeatures();

    my $peptides = $self->getPeptides($peptideType);
    my $maxFeatureNumber = $self->{MaxFeatureNumber};

    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};

	foreach my $pOnePosition (keys %$pOnePositions){
	    
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $featureString = "";
	    my $structureAccessibilityString = $self->getResultField($seqInfo, "Peptide Predicted Accessibility Fraction");
	    my $counter = 1;

	    my @accessibilityFractions = split(',', $structureAccessibilityString);
	    foreach my $accessibilityFraction (@accessibilityFractions){

		my $currentFeatureNumber = $maxFeatureNumber + $counter;
		$featureString .= $currentFeatureNumber . ":$accessibilityFraction ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{structureAccessibilityFeature} = $featureString;
	}
    }
    my $peptideLength = $self->{PeptideLength};
    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}


sub printPeptides{
    my ($self, $peptideType) = @_;

    my $peptides = $self->getPeptides($peptideType);
    my $maxFeatureNumber = $self->{MaxFeatureNumber};

    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};

	foreach my $pOnePosition (keys %$pOnePositions){
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $accession = $self->getResultField($seqInfo, "Uniprot Accession");
#	    $self->writeLog("next string: " . $accession . ":");
	    print STDERR "next string: " . $accession . " seq id $modbaseSeqId p1 $pOnePosition status $peptideType:\n";

	    my $allFeatureStrings = $seqInfo->{allFeatures};
	    foreach my $featureString (keys %$allFeatureStrings){
		my $value = $allFeatureStrings->{$featureString};
#		$self->writeLog("\t$featureString\t$value");
		print STDERR "\t$featureString\t$value\n"; 
	    }
	}
    }
}




sub getResultField{
    
    my ($self, $seqInfo, $resultField) = @_;
    
    my $resultFieldStrings = $self->{ResultFieldStrings};
    if ($resultFieldStrings->{$resultField}){
	
	my $result = $seqInfo->{$resultField};
	return $result;
    }
    else {
	my $allResultFieldStrings = join (", ", keys %$resultFieldStrings);
	print STDERR "Error:  result field $resultField is not a valid protease pipeline result\n";
	print STDERR "Valid result fields include the following:\n";
	print STDERR $allResultFieldStrings . "\n";
	exit(1);
    }
}

sub trainModel{

    my ($self, $trainingSet, $trainingSetFileName, $modelFileName) = @_;

    $self->writeTrainingSetFile($trainingSet, $trainingSetFileName);

    $self->trainSvm($trainingSetFileName, $modelFileName);

}


sub trainSvm{
    my ($self, $trainingSetFileName, $modelFileName) = @_;

    #set RBF params
    my $c = $self->getParam("c");
    my $gamma = $self->getParam("gamma");
    my $flags = "";

    if ($c > 0 && $gamma > 0){
	$flags = "-g $gamma -c $c";
    }

    #get command for running
    my $svmSoftware = $self->getParam('svm_software');
    my $svmCmdName;
    if ($svmSoftware eq 'svmlite'){
	$svmCmdName = $self->getParam("svm_learn_command");
    }
    elsif ($svmSoftware eq 'libsvm'){
#	$svmCmdName = "$pipelineDir/bin/svm-train -b 1";   #only used this in benchmarking for paper, can put it back in if there is a demand for it
    }
    else {
	print STDERR "error: did not get expected svm software name (expected svmlite or libsvm, got $svmSoftware)\n";
	exit(1);
    }

    #train model
    my $svmCmd = "$svmCmdName $flags $trainingSetFileName $modelFileName";
    print STDERR "running svm train command $svmCmd\n";
    
    system($svmCmd);
    
}


sub writeTrainingSetFile{
    my ($self, $trainingSet, $trainingSetFileName) = @_;

    my $trainingSetFh = FileHandle->new(">" . $trainingSetFileName) || die "could not open $trainingSetFileName file for training set writing\n";
    
    #go through all sequences, write one line for each peptide
    foreach my $sequenceId (keys %$trainingSet){
	my $startPositions = $trainingSet->{$sequenceId};
	foreach my $startPosition (keys %$startPositions){
	    my $status = $startPositions->{$startPosition};
	    my $seqInfo =  $self->{Peptides}->{$status}->{$sequenceId}->{$startPosition};

	    #these are feature strings encoded for SVM input (in createPeptideFeatures())
	    my $allFeatureStrings = $seqInfo->{allFeatures};
	    my $outputLine = "";
	    
	    unless ($allFeatureStrings){
		print STDERR "ERROR: Did not find sequence id $sequenceId P1 position $startPosition status $status, retrieved from  training/test set, in data loaded from Protease Pipeline results\n";
		exit(1);  #TODO -- figure out how to handle cluster errors
	    }
	    
	    #prepend line with classification
	    if ($status eq $POSITIVE){
		$outputLine .= "1 ";
	    }
	    if ($status eq $NEGATIVE){
		$outputLine .= "-1 " ;
	    }
		
	    #join all features into one string, write to line
	    my @featureStringArray = sort({$a <=> $b} values %$allFeatureStrings);
	    foreach my $entry (@featureStringArray){
		if ($entry =~ /\S+/){
		    $outputLine .= $entry . " ";
		}
	    }
	    print $trainingSetFh $outputLine . "\n";
	}
    }
}


sub testModel{

    my ($self, $testSet, $testSetFileName, $modelFileName) = @_;

    #test set id list: just a list where each entry is seq_id start_position classification
    my $testSetIdList = $self->writeTestSetFile($testSet, $testSetFileName);

    $self->testSvm($testSetFileName, $modelFileName);
    
    return $testSetIdList;

}

sub writeTestSetFile{
    my ($self, $testSet, $testSetFileName) = @_;
    my $testSetFh = FileHandle->new(">" . $testSetFileName) || die "could not open $testSetFileName file for test set writing\n";
    my @testSetIdList;

    foreach my $sequenceId (keys %$testSet){
	next unless ($sequenceId);  #this is only an issue with homolog map, can probably remove
	my $startPositions = $testSet->{$sequenceId};
	foreach my $startPosition (keys %$startPositions){
	    my $status = $startPositions->{$startPosition};
	    my $seqInfo =  $self->{Peptides}->{$status}->{$sequenceId}->{$startPosition};
	    my $allFeatureStrings = $seqInfo->{allFeatures};
	    my $outputLine = "";
	    unless ($allFeatureStrings){
		print STDERR "ERROR: Did not find sequence id $sequenceId P1 position $startPosition, status $status retrieved from test set, in data loaded from Protease Pipeline results\n";
		exit(1);
	    }
	    my @featureStringArray = sort({$a <=> $b} values %$allFeatureStrings);

	    #prepend line with 0 for unknown
	    $outputLine = "0 ";
		
	    foreach my $entry (@featureStringArray){
		if ($entry =~ /\S+/){
		    $outputLine .= $entry . " ";
		}
	    }
	    print $testSetFh $outputLine . "\n";
	    my $peptideTestSetId = "$sequenceId $startPosition $status";  #save test set ID in list so we know the order they were written (necessary when processing results)
	    push (@testSetIdList, $peptideTestSetId);
	}
    } 
    return \@testSetIdList;
}

sub testSvm{

    my ($self, $testSetFileName, $modelFile) = @_;

    #get svm output file name (list of scores, one per line)
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $runName = $self->getParam("run_name");
    my $svmScoreFileName = $self->getParam("svm_score_file_name");
    my $fullScoreFileName = $pipelineDir . "/" . $runName . "/" . $svmScoreFileName;
    
    #get svm software to use
    my $svmSoftware = $self->getParam('svm_software');
    my $svmCmdName;

    if ($svmSoftware eq 'svmlite'){
	$svmCmdName = $self->getParam("svm_classify_command");
    }
    elsif ($svmSoftware eq 'libsvm'){
	$svmCmdName = "$pipelineDir/bin/svm-predict -b 1";  #only used this in benchmarking for paper, can put it back in if there is a demand for it
    }
   
    #run command
    my $applySvmCmd = "$svmCmdName $testSetFileName $modelFile $fullScoreFileName";
    print STDERR "svm model: applying svm to application set with command $applySvmCmd\n";
    system($applySvmCmd);
}

sub setNoPeptidesParsed{
    my ($self) = @_;
    $self->{noPeptidesParsed} = 1;
}
sub getNoPeptidesParsed{
    my ($self) = @_;
    return $self->{noPeptidesParsed};
}

sub parseResults{

    my ($self, $testSetIdList) = @_;

    #create $testSetLines which is keyed on index of the peptide in the order as it was written to the test set file
    my $testSetLines;
    my $testSetLineCounter = 0;
    foreach my $testSetEntry (@$testSetIdList){
	$testSetLineCounter++;
	$testSetLines->{$testSetLineCounter} = $testSetEntry;
    }
 
    #get score file name
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $runName = $self->getParam("run_name");
    my $svmScoreFileName = $self->getParam("svm_score_file_name");
    my $fullScoreFileName = $pipelineDir . "/" . $runName . "/" . $svmScoreFileName;
   
    my $scoreFh = FileHandle->new("<" . $fullScoreFileName) || die "could not open $fullScoreFileName\n";
    my $scoreLines;
    my $scoreLineCounter = 0;
   
    my $svmSoftware = $self->getParam('svm_software');
    
    #get score in each line, save in $scoreLines which is keyed on line number
    while (<$scoreFh>){
	chomp;
	my $line = $_;

	if ($svmSoftware eq 'libsvm'){
	    next if ($line =~ /label/);
	    my @cols = split('\s', $line);
	    $line = $cols[2];
	}
	$scoreLineCounter++;
	$scoreLines->{$scoreLineCounter} = $line;
    }

    my $results;

    #merge the two hashes into $results, which is of form $results->{$peptideId}->{score} = $score
    foreach my $lineNumber (sort ({$scoreLines->{$a} <=> $scoreLines->{$b}}  keys %$scoreLines)){   #TODO -- figure out why we are sorting these if we are just putting it into hash
	my $score = $scoreLines->{$lineNumber};
	my $input = $testSetLines->{$lineNumber};
	my $id;
	if ($input =~ /(\S+\s\S+\s\S+)/){  #id is of the form seq_id start_position classification TODO -- figure out why we are regexping it, should be the only thing in the entry (maybe needed for other SVM software?)
	    $id = $1;
	}
	else {
	    die "did not parse proper input from SVM result line $input\n";
	}
	$results->{$id}->{score} = $score;
    
    }  
    return $results;
}


sub addResultsToPeptideInfo{

    my ($self, $resultSet, $type) = @_;
    
    my $peptides = $self->getPeptides($type);
    foreach my $peptideSequenceId (keys %$peptides){
	my $pOnePositions = $peptides->{$peptideSequenceId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $idLine = $peptideSequenceId . " " . $pOnePosition . " " . $type;   #use unique identifiers to recreate identification string that was output in SVM results file
	    my $score = $resultSet->{$idLine}->{score};	 
	    $self->{Peptides}->{$type}->{$peptideSequenceId}->{$pOnePosition}->{"SVM Score"} = $score;
	    
	    unless ($score) {
		print STDERR "error: did not get score for id line $idLine\n";
	    }	
	}
    }
}

sub addBenchmarkScoresToPeptideInfo{

    my ($self, $type) = @_;
    my $peptides = $self->getPeptides($type);
    print STDERR "getting benchmark rates for all scores\n";
    my $ratesAtEachScore = $self->getBenchmarkRatesForAllScores();
    foreach my $peptideSequenceId (keys %$peptides){
	my $pOnePositions = $peptides->{$peptideSequenceId};
	foreach my $pOnePosition (keys %$pOnePositions){
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $score = $seqInfo->{"SVM Score"};
	    
	    my ($fpr, $tpr) = $self->getBenchmarkRatesForOneScore($ratesAtEachScore, $score);

	    $self->{Peptides}->{$type}->{$peptideSequenceId}->{$pOnePosition}->{"TPR"} = $tpr;
	    $self->{Peptides}->{$type}->{$peptideSequenceId}->{$pOnePosition}->{"FPR"} = $fpr;
	}
    }
}

sub getBenchmarkRatesForAllScores{
    my ($self) = @_;

    my $benchmarkScoreFile = $self->getParam("benchmark_score_file");
    print STDERR "reading benchmark score file $benchmarkScoreFile for results\n";
    my $benchmarkScoreFh = FileHandle->new("<" . $benchmarkScoreFile) || die "could not open benchmark score file  $benchmarkScoreFile\n";
    my $scoreBins;
    
    #read benchmark file (this is just copied from the results of the SVM run)
    while (<$benchmarkScoreFh>){
	chomp;
	my $line = $_;
	my ($fpRate, $tpRate, $svmScore) = split('\t', $line);
	$svmScore = sprintf('%.3f', $svmScore);
	$scoreBins->{$svmScore}->{tpRate} = $tpRate;
	$scoreBins->{$svmScore}->{fpRate} = $fpRate;
    }
    return $scoreBins;
}


sub getBenchmarkRatesForOneScore{
    my ($self, $ratesAtEachScore, $svmScore) = @_;
 
    my @allScores = sort ({$a <=> $b} keys %$ratesAtEachScore);  #sorts from worst SVM score (most negative) to best
    my $worstScore = $allScores[0];
    $svmScore = sprintf('%.3f', $svmScore);
    #get the tp rate and fp rate at this score.  Since not all scores will be represented, if we don't see one, 
    #just go down to the first one we do see (will only be the case for the extreme score values)
    
    my $results = $ratesAtEachScore->{$svmScore};
    while (!($results)){     
	$svmScore -= .001;
	$svmScore = sprintf('%.3f', $svmScore);   

	$results = $ratesAtEachScore->{$svmScore};
	if ($svmScore < $worstScore){ #we got something so bad that it is below the benchmark, just take the worst thing 
	    $results = $ratesAtEachScore->{$worstScore};
	}
    }
    unless ($results){
	die "did not get any benchmark results for svm score $svmScore\n";   #might see this because method hasn't been fully testedy
    }
    my $tpRate = $results->{tpRate};
    my $fpRate = $results->{fpRate};
    
    return ($fpRate, $tpRate);
}

sub getPeptides{
    my ($self, $peptideType) = @_;
    my $peptides = $self->{Peptides}->{$peptideType};
    return $peptides;
}


sub initialize{
    my ($self) = @_;

    $self->makeResidueCodes();
    $self->{MaxFeatureNumber} = 0;
    my $peptideLength = $self->getParam("peptide_length");
    $self->{PeptideLength} = $peptideLength;
}


sub makeResidueCodes{
    
    my ($self) = @_;
    
    my @residueList = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");
    #each residue is represented by an array of size 20 where every entry is 0 except for entry indexed by the position of the residue in @residueList, which is set to 1
    #implicit/explicit zeroes are handled when sequence features are created
    my $currentResidueCount = 0;
    foreach my $residue (@residueList){

	my @array = (0);
	for (my $i = 1; $i < 20; $i++){
	    push(@array, 0);
	}
	$array[$currentResidueCount] = 1;
	$self->{ResidueStrings}->{$residue} = \@array;

	$currentResidueCount++;
    }
}

sub useStructureFeatures{
    my ($self) = @_;
    my $useStructureFeatures = $self->getParam("use_structure_features");
    if ($useStructureFeatures eq "no"){
	return 0;
    }
    elsif ($useStructureFeatures eq "yes"){
	return 1;
    }
    else {
	print STDERR "ERROR: please set parameter 'use_structure_features' to either yes or no\n";
	exit(1);
    }

}
sub trainsOnNegatives{
    return 1;
}

sub resetFeatureNumber{

    my ($self) = @_;

    $self->{MaxFeatureNumber} = 0;
}


return 1;
