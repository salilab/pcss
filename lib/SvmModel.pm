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

############################################################################################################################################################     
# loadPeptides
# 
# Reads peptidePipelineResults.txt file output by PeptidePipeline.pm run that evaluated all structural features for the input peptides. Loads all features
# into data structures.
#
# Peptides are all stored in $self->{Peptides}.  This is a large hash of the following format:
# $self->{Peptides}->{$peptideType}->{$modbaseSeqId}->{$peptideStartPosition} = $seqInfo
#
# $peptideType is the Classification of the peptide ("positive", "negative", "Application"); this is set by the frontend and stored in peptide result file.
# $seqInfo is another hash where the keys are the column names as read from the peptide result file and the values are the values for that column for the peptide.
#
# The Model also saves the number of peptides it loaded of each type by calling $self->setPeptideCount().
#
# PARAM  $peptideFileName: the name of the result file from which to read all peptide info.
# RETURN $sequenceCounts: a hash containing unique modbase sequence ids loaded from the results file, of the form
#        $sequenceCounts->{$peptideType}->{$modbaseSeqId} = 1
# RETURN $peptideCoutns: a hash where the keys are the different peptide classifications and the values are the number of peptides loaded with that classification.
############################################################################################################################################################     
sub loadPeptides{

    my ($self, $peptideFileName) = @_;

    my $peptideCounts;    
    my $sequenceCounts;
    my $columnMap;    

    my $peptideFh = FileHandle->new("<" . $peptideFileName);
    unless ($peptideFh =~ /FileHandle/){
	$self->writeError("file_missing", "loadPeptides", "could not open PeptidePipeline result file $peptideFileName: $!");
	return 0;
    }
    my $firstLine = <$peptideFh>;
   
    while ($firstLine =~ /^\s+$/){
	$firstLine = <$peptideFh>;
    }

    #if PeptidePipeline resulted in output error (e.g. column names couldn't be written), it is specified on the first line
    my $keywordOutputError = $self->getParam("keyword_output_error");
    if ($firstLine =~ /$keywordOutputError/){
	$self->writeError($keywordOutputError, "loadPeptides", "Read in output error keyword ($keywordOutputError) when reading in PeptidePipeline feature result file");
	return 0;
    }
    else {
	#Read in column names from first line; make map of column display order to column name
	my @cols = split('\t', $firstLine);
	my $i = 0;
	foreach my $col (@cols){
	    $self->{ResultFieldStrings}->{$col} = 1;
	    $columnMap->{$i} = $col;
	    $i++;
	}
    }

    #read rest of result file
    my $onFirstContentLine = 0;
    while (<$peptideFh>){
		
	chomp;
	my $line = $_;

	#first content line is the first one below the column headers. Other global errors are written here
	if ($onFirstContentLine == 0){
	    $onFirstContentLine = 1;
	    my $globalErrorRegex = $self->getGlobalErrorRegex();

	    #check if other errors occured in previous PeptidePipeline feature run
	    if ($line =~ /($globalErrorRegex)/){
		my $errorName = $1;
		$self->writeError($errorName, "loadPeptides", "Read in global error keyword ($errorName) when reading in PeptidePipeline feature result file");
		return 0;
	    }
    	}
	
	#process result line with all features for peptide; add to $self->{Peptides}
	my @cols = split('\t', $line);
	
	my $currentSeqInfo;
	my $i = 0;
	foreach my $colValue (@cols){     #read in field values 
	    my $currentColumn = $columnMap->{$i};
	    $currentSeqInfo->{$currentColumn} = $colValue;
	    $i++;
	}
	
	#check if no peptides were parsed for this sequence in the previous PeptidePipeline feature run
	my $noPeptidesParsedKeyword = $self->getParam("keyword_no_peptides_parsed");
	if ($line =~ /$noPeptidesParsedKeyword/){
	    my $seqId = $currentSeqInfo->{"Sequence ID"};
	    my $uniprotAccession = $currentSeqInfo->{"Uniprot Accession"};
	    $self->writeNoPeptidesParsedForSeq($seqId, $uniprotAccession);
	}
	else {
	    #get relevant info and add to $self->{Peptides}
	    my $seqId = $currentSeqInfo->{"Sequence ID"};
	    my $startPosition = $currentSeqInfo->{"Peptide Start"};
	    my $peptideType = $currentSeqInfo->{"Classification"};  #TODO -- parameterize this (was set in frontend and is assumed to match one of the global types)
	    $peptideCounts->{$peptideType}++;
	    $sequenceCounts->{$peptideType}->{$seqId} = 1;
	    $self->{Peptides}->{$peptideType}->{$seqId}->{$startPosition} = $currentSeqInfo;   
	}
    }
    
    foreach my $peptideType (keys %$peptideCounts){
	my $peptideCount = $peptideCounts->{$peptideType};
	$self->setPeptideCount($peptideCount, $peptideType);
    }
    return ($sequenceCounts, $peptideCounts);
}




############################################################################################################################################################     
# createPeptideFeatures
#
# Go through all peptides loaded and read sequence and structure info. Convert this into code expected by the SVM software. All features are stored in the 
# hash keyed on $self->{Peptides}->{$peptideType}->{$modbaseSeqId}->{$peptideStartPosition}->{allFeatures} (one for each peptide)
#
# There can be hundreds of features for one peptide. Each feature is represented by its feature number (1 through n, n = total number of features) and then
# its value. The result will be a long string of feature number/value pairs, where the number and value are separated by ':'
#
# A running call to $self->{MaxFeatureNumber}, which is updated after each method run here, lets the Model know what the current feature number is. 
# $self->resetFeatureNumber MUST be called when loading multiple peptide feature datasets (eg, loading first positives and then negatives). 
# (This is not really enforced though as the createXYZfeatures() methods are never called other than in this method).
#
# PARAM  $peptideType: specifies the classification of peptides for which features will be created ("positive", "negative", "Application").
# RETURN NULL
############################################################################################################################################################     
sub createPeptideFeatures{
    my ($self, $peptideType) = @_;

    $self->resetFeatureNumber(); 

    $self->createPeptideSequenceFeatures($peptideType);
    
    $self->createDisorderFeatures($peptideType);

    $self->createStructureAccessibilityFeatures($peptideType);
    
    $self->createSequenceSecondaryStructureFeatures($peptideType);
    
    $self->createStructureSecondaryStructureFeatures($peptideType);
}


############################################################################################################################################################     
# createPeptideSequenceFeatures
# 
# For all peptides of type $peptideType, convert peptide residue sequence into a code readable by the SVM software. Each residue gets a 
# feature number corresponding to its identity and position in the peptide sequence of (20n + i) where n is the 0-based position and i is the alphabetical
# position of the one-letter residue code. The value for a residue is 1. Thus, Asp in the fourth position would be 20(3) + 3 = 63, and the feature would be
# represented as 63:1. If the "zeroes" parameter is set to explicit, then all other residues get their own feature number which is set to 0 (i.e., if Asp is
# in the fourth position, then 60 - 62 and 64-79 are all set to zero explicitly); otherwise, these feature numbers aren't listed.
#
# The full feature string is stored in $seqInfo->{allFeatures}->{peptideSeqFeature} (one $seqInfo for each peptide).
#
# PARAM  $peptideType: specifies the classification of peptides for which features will be created ("positive", "negative", "Application").
# RETURN NULL
############################################################################################################################################################     
sub createPeptideSequenceFeatures{
    
    my ($self, $peptideType) = @_;

    my $peptides = $self->getPeptides($peptideType);
    
    my $maxFeatureNumber = $self->{MaxFeatureNumber};

    my $zeroes = $self->getParam("zeroes");

    my $peptideLength = $self->{PeptideLength};

    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};

	foreach my $pOnePosition (keys %$pOnePositions){

	    my $featureString = "";
	    
	    #get peptide sequence
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $peptideSequence = $self->getResultField($seqInfo, "Peptide Sequence");
	    
	    my @seqArray = split('', $peptideSequence);
    
	    my $counter = 1;
	    my $featureNumber;
	    my $residueCounter = 0;		    
	    foreach my $nextResidue (@seqArray){
		last if ($residueCounter >= $peptideLength);
		$residueCounter++;
		
		#get the @residueArray for this residue, 20 elements where each is 0 except for the element set to 1, indexed by the alphabetical position of $nextResidue.
		my $residueArray = $self->{ResidueStrings}->{$nextResidue};
		
		unless ($residueArray){
		    my $msg = "Did not get a residue array for residue $nextResidue in modbase sequence $modbaseSeqId start position $pOnePosition peptide $peptideSequence ";
		    $msg .= " (likely $nextResidue is not one of the 20 standard amino acids)";
		    $self->writeError("invalid_residue", "createPeptideSequenceFeatures", $msg);
		    return 0;
		}

		foreach my $nextFeature(@$residueArray){
		    $featureNumber = $counter + $maxFeatureNumber;
		    $counter++;
		    
		    #if zeroes not "explicit" then don't write them out
		    next if ($zeroes eq "implicit" && $nextFeature == 0);
		    $featureString .= $featureNumber . ":" . $nextFeature . " ";
		}
	    }
	    $seqInfo->{allFeatures}->{peptideSeqFeature} = $featureString;
	}
    }

    $maxFeatureNumber += $peptideLength * 20;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}


############################################################################################################################################################     
# createDisorderFeatures
#
# The values for disorder features are just the raw Disopred scores. One feature is written for each position in the peptide, with the value 
# corresponding to the disopred score at that position.
#
# PARAM  $peptideType: specifies the classification of peptides for which features will be created ("positive", "negative", "Application").
# RETURN NULL
############################################################################################################################################################     
sub createDisorderFeatures{

    my ($self, $peptideType) = @_;
    return unless $self->useStructureFeatures();

    my $peptideLength = $self->{PeptideLength};
    my $peptides = $self->getPeptides($peptideType);
    my $maxFeatureNumber = $self->{MaxFeatureNumber};
    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};

	foreach my $pOnePosition (keys %$pOnePositions){
	    my $counter = 1;
	    my $featureString = "";
	    my $seqInfo = $pOnePositions->{$pOnePosition};

	    #get Disopred scores for each position, convert to feature string
	    my $disorderString = $self->getResultField($seqInfo, "Disopred Scores");
	    my @disorderScores = split(',', $disorderString);
	    foreach my $disorderScore (@disorderScores){
		last if ($counter > $peptideLength);
		my $currentFeatureNumber = $maxFeatureNumber + $counter;
		$featureString .= $currentFeatureNumber . ":$disorderScore ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{disorderFeature} = $featureString;
	}
    }

    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}


############################################################################################################################################################     
# createDisorderFeatures
#
# The values for sequence-based secondary structure features are just the raw Psipred scores. One feature is written for each position in the peptide, with the value 
# corresponding to the psipred score at that position.
#
# PARAM  $peptideType: specifies the classification of peptides for which features will be created ("positive", "negative", "Application").
# RETURN NULL
############################################################################################################################################################     
sub createSequenceSecondaryStructureFeatures{

    my ($self, $peptideType) = @_;

    return unless $self->useStructureFeatures();

    my $peptides = $self->getPeptides($peptideType);
    
    my $maxFeatureNumber = $self->{MaxFeatureNumber};
    my $peptideLength = $self->{PeptideLength};
    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};
	
	foreach my $pOnePosition (keys %$pOnePositions){
	    my $counter = 1;
	    my $featureString = "";
	    my $seqInfo = $pOnePositions->{$pOnePosition};

	    #get the Psipred score at each position, set as value for the feature
	    my $secondaryStructureString = $self->getResultField($seqInfo, "PSIPRED Scores");
	    my @psipredScores = split(',', $secondaryStructureString);

	    foreach my $psipredScore (@psipredScores){
		last if ($counter > $peptideLength);
		my $currentFeatureNumber = $maxFeatureNumber + $counter;  
		$featureString .= $currentFeatureNumber . ":$psipredScore ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{sequenceSecondaryStructureFeature} = $featureString;
	}
    }

    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}

############################################################################################################################################################     
# createStructureSecondaryStructureFeatures
#
# The values for structure-based secondary structure features were assigned in the PeptidePipeline step. They are simply 1,2,3, depending on structure type
# (loop, helix, sheet) assigned by DSSP. One feature is written for each position in the peptide, with the value corresponding to the (1,2,3) assignment at 
# that position
#
# PARAM  $peptideType: specifies the classification of peptides for which features will be created ("positive", "negative", "Application").
# RETURN NULL
############################################################################################################################################################     
sub createStructureSecondaryStructureFeatures{

    my ($self, $peptideType) = @_;
    return unless $self->useStructureFeatures();

    my $peptides = $self->getPeptides($peptideType);
    
    my $maxFeatureNumber = $self->{MaxFeatureNumber};
    my $peptideLength = $self->{PeptideLength};
    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};
 
	foreach my $pOnePosition (keys %$pOnePositions){
	    
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $featureString = "";

	    #get structure calls and assign as features
	    my $structureSecondaryStructureString = $self->getResultField($seqInfo, "Peptide Structure Values");
	    my $counter = 1;

	    my @secondaryStructureValues = split(',', $structureSecondaryStructureString);
	    foreach my $secondaryStructureValue (@secondaryStructureValues){
		last if ($counter > $peptideLength);
		my $currentFeatureNumber = $maxFeatureNumber + $counter;
		$featureString .= $currentFeatureNumber . ":$secondaryStructureValue ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{structureSecondaryStructureFeature} = $featureString;
	}
    }

    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}


############################################################################################################################################################     
# createStructureAccessibilityFeatures
#
# The values for structure-based solvent accessibilty features are just the accessibility fraction as assessed by DSSP. One feature is written for each 
# position in the peptide, with the value corresponding to the fraction assignment at that position
#
# PARAM  $peptideType: specifies the classification of peptides for which features will be created ("positive", "negative", "Application").
# RETURN NULL
############################################################################################################################################################     
sub createStructureAccessibilityFeatures{

    my ($self, $peptideType) = @_;
    return unless $self->useStructureFeatures();

    my $peptides = $self->getPeptides($peptideType);
    my $maxFeatureNumber = $self->{MaxFeatureNumber};
    my $peptideLength = $self->{PeptideLength};
    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};

	foreach my $pOnePosition (keys %$pOnePositions){
	    
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $featureString = "";

	    #get accessibility fraction and assign as value for the next feature
	    my $structureAccessibilityString = $self->getResultField($seqInfo, "Peptide Predicted Accessibility Fraction");
	    my $counter = 1;

	    my @accessibilityFractions = split(',', $structureAccessibilityString);
	    foreach my $accessibilityFraction (@accessibilityFractions){
		last if ($counter > $peptideLength);
		my $currentFeatureNumber = $maxFeatureNumber + $counter;
		$featureString .= $currentFeatureNumber . ":$accessibilityFraction ";
		$counter++;
	    }
	    $seqInfo->{allFeatures}->{structureAccessibilityFeature} = $featureString;
	}
    }

    $maxFeatureNumber += $peptideLength;
    $self->{MaxFeatureNumber} = $maxFeatureNumber;
}


############################################################################################################################################################     
# printPeptides
# 
# prints all peptides with the classification specified by $peptideType. Simple output / debug function; output could probably be cleaned up.
# Peptides are printed to the log file.
############################################################################################################################################################     
sub printPeptides{
    my ($self, $peptideType) = @_;

    my $peptides = $self->getPeptides($peptideType);

    foreach my $modbaseSeqId (keys %$peptides){
	my $pOnePositions = $peptides->{$modbaseSeqId};

	foreach my $pOnePosition (keys %$pOnePositions){
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $accession = $self->getResultField($seqInfo, "Uniprot Accession");
	    $self->writeLog("next string: " . $accession . " seq id $modbaseSeqId p1 $pOnePosition status $peptideType");

	    my $allFeatureStrings = $seqInfo->{allFeatures};
	    foreach my $featureString (keys %$allFeatureStrings){
		my $value = $allFeatureStrings->{$featureString};
		$self->writeLog("\t$featureString\t$value");
	    }
	}
    }
}


############################################################################################################################################################     
# getResultField
#
# Gets the value specified by $resultField from the hash $seqInfo. Just a simple hash lookup, but since these are all result fields specified by columns
# and other output files, provides some protection against looking for an expected value in $seqInfo using a random string.
#
# PARAM  $seqInfo: Hash specifying all peptide features; keys are column display names read in from the PeptidePipeline results file (or otherwise set in this
#                  pipeline) and values are the values for that feature for the peptide.
# PARAM  $resultField: the name of the feature (also should be a column name) to look for in $seqInfo.
# RETURN value at $seqInfo->{$resultField} if it exists; otherwise throws an error
############################################################################################################################################################     
sub getResultField{
    
    my ($self, $seqInfo, $resultField) = @_;
    
    my $resultFieldStrings = $self->{ResultFieldStrings};
    if ($resultFieldStrings->{$resultField}){
	
	my $result = $seqInfo->{$resultField};
	return $result;
    }
    else {
	my $allResultFieldStrings = join (", ", keys %$resultFieldStrings);
	my $errorMessage =  "Error:  result field $resultField is not a valid protease pipeline result. Valid result fields include the following:";
	$errorMessage .= $allResultFieldStrings;
	$self->writeError("invalid_result_field", "getResultField", $errorMessage); 
	#called in too many places to try to return an error and handle in all cases. However, $self->writeError() calls setError() which should
	#kill pipeline processing in future steps.
    }
}

############################################################################################################################################################     
# trainModel()
# 
# Writes peptide features to an SVM-readable file, and then creates an actual model using the SVM software.
# 
# RETURN NULL
############################################################################################################################################################     

sub trainModel{

    my ($self, $trainingSet, $trainingSetFileName, $modelFileName) = @_;

    my $result = $self->writeTrainingSetFile($trainingSet, $trainingSetFileName);
    return 0 unless $result;

    $self->trainSvm($trainingSetFileName, $modelFileName);

}

############################################################################################################################################################
# trainSvm
#
# Prepare command to run the SVM training executable and then run it.
# 
# PARAM  $trainingSetFileName: name of the file that will be used as input to SVM software
# PARAM  $modelFileName: name of the SVM model (also input to SVM software) to output and use later
# RETURN NULL
############################################################################################################################################################
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
	#$svmCmdName = "$pipelineDir/bin/svm-train -b 1";   #only used this in benchmarking for paper, can put it back in if there is a demand for it
    }

    #train model
    my $svmCmd = "$svmCmdName $flags $trainingSetFileName $modelFileName";
    $self->writeLog("training svm with command $svmCmd");
    my $output = `$svmCmd 2>&1`;

    my $exitCode = $?;
    unless ($exitCode == 0){
	$self->writeError("svm_failure", "trainSvm", "SVM train command $svmCmd exited with non-zero exit code ($exitCode)");
	$self->writeLog("output of SVM cmd: $output");
    }
}


############################################################################################################################################################     
# writeTrainingSetFile()
# Write out peptide feature file that will be read by SVM software. Each peptide gets its own line that is appropriately coded as feature number / value list
# (see createPeptideFeatures() for description). Lines are preceeded with either "1" or "-1" to indicate that their classification
# PARAM  $trainingSet: hash specifying the sequences and peptides that the SVM will score. 
#                      format: $testSet->{$modbaseSequenceId}->{$peptideStartPosition} = $TYPE ($POSITIVE or $NEGATIVE)
# PARAM  $trainingSetFileName: name of the file that will be used as input to SVM software
# RETURN 1 if success; 0 if errors occur
############################################################################################################################################################     
sub writeTrainingSetFile{
    my ($self, $trainingSet, $trainingSetFileName) = @_;

    my $trainingSetFh = FileHandle->new(">" . $trainingSetFileName);
    unless ($trainingSetFh =~ /FileHandle/){
	$self->writeError("file_missing", "writeTrainingSetFile", "could not open model training input file $trainingSetFileName for writing: $!");
	return 0;
    }
    
    #go through all sequences, write one line for each peptide
    foreach my $sequenceId (keys %$trainingSet){
	my $startPositions = $trainingSet->{$sequenceId};
	foreach my $startPosition (keys %$startPositions){
	    my $status = $startPositions->{$startPosition};
	    my $seqInfo =  $self->{Peptides}->{$status}->{$sequenceId}->{$startPosition};

	    #these are feature strings encoded for SVM input (in createPeptideFeatures())
	    my $allFeatureStrings = $seqInfo->{allFeatures};
	    my $outputLine = "";
	    
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

    return 1;
}


############################################################################################################################################################     
# testModel
# 
# Write out peptide feature file that will be read by SVM software, and then call system command to apply the SVM to score the peptides in the file
# 
# PARAM  $testSet: hash specifying the sequences and peptides that the SVM will score. 
#                  format: $testSet->{$modbaseSequenceId}->{$peptideStartPosition} = $PROTEOME
# PARAM  $testSetFileName: name of the file that will be used as input to SVM software
# PARAM  $modelFileName: name of the SVM model (also input to SVM software) that was trained and benchmarked by the server in a separate training run
# RETURN $testSetIdList: array ref where each entry is of the format "$modbaseSequenceId $startPosition $status". Used to track order in which the peptides
#                        are written to $testSetFileName (necessary when reading results file later)
############################################################################################################################################################     
sub testModel{  #TODO -- change to applyModel

    my ($self, $testSet, $testSetFileName, $modelFileName) = @_;

    my $testSetIdList = $self->writeTestSetFile($testSet, $testSetFileName);

    return 0 unless $testSetIdList;

    $self->testSvm($testSetFileName, $modelFileName);
    
    return $testSetIdList;
}


############################################################################################################################################################ 
# writeTestSetFile
# Write out peptide feature file that will be read by SVM software. Each peptide gets its own line that is appropriately coded as feature number / value list
# (see createPeptideFeatures() for description). Lines are preceeded with "0" to indicate that this is an unclassified peptide.
#
# PARAM  $testSet: hash specifying the sequences and peptides that the SVM will score. 
#                  format: $testSet->{$modbaseSequenceId}->{$peptideStartPosition} = $PROTEOME
# PARAM  $testSetFileName: name of the file that will be used as input to SVM software
# RETURN $testSetIdList: array ref where each entry is of the format "$modbaseSequenceId $startPosition $status". Used to track order in which the peptides
#                        are written to $testSetFileName (necessary when reading results file later)
############################################################################################################################################################
sub writeTestSetFile{
    my ($self, $testSet, $testSetFileName) = @_;
    my $testSetFh = FileHandle->new(">" . $testSetFileName);
    unless ($testSetFh =~ /FileHandle/){
	$self->writeError("file_missing", "writeTestSetFile", "could not open model application input file $testSetFileName for writing: $!");
	return 0;
    }
    my @testSetIdList;

    $self->writeLog("writing peptides to SVM application file $testSetFileName");

    foreach my $sequenceId (keys %$testSet){
	next unless ($sequenceId);  #TODO - this is only an issue with homolog map, can probably remove
	my $startPositions = $testSet->{$sequenceId};
	foreach my $startPosition (keys %$startPositions){

	    #get features for this peptide
	    my $status = $startPositions->{$startPosition};
	    my $seqInfo =  $self->{Peptides}->{$status}->{$sequenceId}->{$startPosition};  
	    my $allFeatureStrings = $seqInfo->{allFeatures};
	    
	    my @featureStringArray = sort({$a <=> $b} values %$allFeatureStrings);

	    #create one line for each peptide. Prepend line with 0 for unknown
	    my $outputLine = "0 ";
		
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


############################################################################################################################################################
# testSvm
#
# Prepare command to run the SVM application executable and then run it.
# 
# PARAM  $testSetFileName: name of the file that will be used as input to SVM software
# PARAM  $modelFileName: name of the SVM model (also input to SVM software) that was trained and benchmarked by the server in a separate training run
# RETURN NULL
############################################################################################################################################################
sub testSvm{

    my ($self, $testSetFileName, $modelFile) = @_;

    #get svm output file name (list of scores, one per line)
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $runName = $self->getParam("run_name");
    my $svmScoreFileName = $self->getParam("svm_score_file_name");
    my $fullScoreFileName = $pipelineDir . "/" . $runName . "/" . $svmScoreFileName;

    unless (-e $modelFile){
	$self->writeError("file_missing", "testSvm", "File $modelFile does not exist");
	return 0;
    }
    
    #get svm software to use
    my $svmSoftware = $self->getParam('svm_software');
    my $svmCmdName;

    if ($svmSoftware eq 'svmlite'){
	$svmCmdName = $self->getParam("svm_classify_command");
    }
    elsif ($svmSoftware eq 'libsvm'){A
	$svmCmdName = "$pipelineDir/bin/svm-predict -b 1";  #only used this in benchmarking for paper, can put it back in if there is a demand for it
    }
   
    #run command
    my $applySvmCmd = "$svmCmdName $testSetFileName $modelFile $fullScoreFileName";
    $self->writeLog("applying svm to application set with command $applySvmCmd");
    my $output = `$applySvmCmd 2>&1`;

    my $exitCode = $?;
    unless ($exitCode == 0){
	$self->writeError("svm_failure", "testSvm", "SVM application command $applySvmCmd exited with non-zero exit code ($exitCode)");
	$self->writeLog("output of SVM cmd: $output");
    }
    
}


############################################################################################################################################################
# parseResults
# 
# Read the results of the SVM application step. Reads the SVM output file which is one line per peptide, where the only thing on the line is the score.
# These are written in the same order as they are listed in $testSetIdList. Use that to map the peptide to the score and return a data structure containing
# all results.
#
# PARAM  $testSetIdList: array ref where each entry is of the format "$modbaseSequenceId $startPosition $status". 
# RETURN $results: hash mapping peptides to scores. This is of the format
#                  $results->{$id}->{score} = $score.
#                  Where $id is the same as $testSetIdList entries, and $score is read from the file.
############################################################################################################################################################
sub parseResults{

    my ($self, $testSetIdList) = @_;

    #create $testSetLines hash, which is keyed on index of the peptide in the order as it was written to the test set file (starting with 1)
    my $testSetLines;
    my $testSetLineCounter = 0;
    foreach my $testSetEntry (@$testSetIdList){
	$testSetLineCounter++;
	$testSetLines->{$testSetLineCounter} = $testSetEntry;
    }
 
    #get score file name. This file has one line written per peptide, which is ONLY the score
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $runName = $self->getParam("run_name");
    my $svmScoreFileName = $self->getParam("svm_score_file_name");
    my $fullScoreFileName = $pipelineDir . "/" . $runName . "/" . $svmScoreFileName;
    $self->writeLog("reading $fullScoreFileName");
    my $scoreFh = FileHandle->new("<" . $fullScoreFileName);
    unless ($scoreFh =~ /FileHandle/){
	$self->writeError("file_missing", "parseResults", "could not open svm score file $fullScoreFileName: $! Check the logs for possible errors in calling the svm executable");
	return 0;  #this is tested, but not in the regression test since I couldn't figure out how to break it (file is dynamically generated). Simple enough though and most
	           #of the time an error that results in no score file will be caught due to the svm software not exiting correctly, including if it couldn't write to the score file.
    }

    my $scoreLines;
    my $scoreLineCounter = 0;
   
    my $svmSoftware = $self->getParam('svm_software');
    
    #get score from each line, save in $scoreLines which is a hash of line number to score
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

    #merge the two hashes of line number to peptide, and line number to score, into $results, which is of form $results->{$peptideId}->{score} = $score
    foreach my $lineNumber (sort ({$scoreLines->{$a} <=> $scoreLines->{$b}}  keys %$scoreLines)){   
	my $score = $scoreLines->{$lineNumber};
	my $input = $testSetLines->{$lineNumber};
	my $id;
	if ($input =~ /(\S+\s\S+\s\S+)/){  
	    $id = $1;
	}
	else {
	    
	    $self->writeError("results_error", "parseResults", "Did not get expected line format for peptide id (expected 'seq_id start_position classification'; got $input");
	    return 0; #tested this, but not in regression test (couldn't figure out how to break it)
	}
	$results->{$id}->{score} = $score;
    }  
    return $results;
}

############################################################################################################################################################
# addResultsToPeptideInfo
# 
# Given peptides held by the model, and results of SVM application, add the score to the peptides $seqInfo hash containing all the other features.
# ($seqInfo located in $self->{Peptides}->{$peptideType}->{$modbaseSeqId}->{$peptideStartPosition} = $seqInfo)
#
# PARAM  $resultSet: hash returned by $self->parseResults(), of the format $results->{$id}->{score} = $score.
#                    where $id is of the format  "$modbaseSequenceId $startPosition $status"
# PARAM  $type:      classification of the peptides to which we are adding scores ("positive", "negative", "application")
# RETURN 1 if normal flow, 0 if the score for a peptide was not found
############################################################################################################################################################
sub addResultsToPeptideInfo{

    my ($self, $resultSet, $type) = @_;
    my $peptides = $self->getPeptides($type);
    foreach my $peptideSequenceId (keys %$peptides){
	my $startPositions = $peptides->{$peptideSequenceId};
	
	foreach my $startPosition (keys %$startPositions){
	    my $seqInfo = $startPositions->{$startPosition};
	    my $idLine = $peptideSequenceId . " " . $startPosition . " " . $type;   #use unique identifiers to recreate identification string that was output in SVM results file
	    my $score = $resultSet->{$idLine}->{score};	 
	    $self->{Peptides}->{$type}->{$peptideSequenceId}->{$startPosition}->{"SVM Score"} = $score;
	    
	    unless ($score) {
		my $msg = "Did not get score for sequence $peptideSequenceId start position $startPosition using id line $idLine";
		$msg .= "Check the svm score output file to make sure there is an entry for every peptide and that none of the scores is 0";
		$self->writeError("results_error", "addResultsToPeptideInfo", $msg);
		return 0;
	    }  #tested, but not in regression test	
	}
    }
    return 1;
}

############################################################################################################################################################
# addBenchmarkScoresToPeptideInfo
# 
# For each peptide, read benchmark result file that was output in training mode and see where the newly generated SVM scores would fall in that benchmark file,
# in terms of the false positive and true positive rates they would have resulted in. Save in "TPR" and "FPR" feature values in the peptide's $seqInfo hash.
# 
# PARAM  $type: classification of the peptides to process ('positive', 'negative', 'application')
# RETURN NULL
############################################################################################################################################################
sub addBenchmarkScoresToPeptideInfo{

    my ($self, $type) = @_;
    my $peptides = $self->getPeptides($type);

    #load all benchmark scores from file
    my $ratesAtEachScore = $self->getBenchmarkRatesForAllScores();
    return 0 unless $ratesAtEachScore;

    foreach my $peptideSequenceId (keys %$peptides){
	my $pOnePositions = $peptides->{$peptideSequenceId};
	foreach my $pOnePosition (keys %$pOnePositions){
	    my $seqInfo = $pOnePositions->{$pOnePosition};
	    my $score = $seqInfo->{"SVM Score"};
	    
	    #see where this score falls in overall benchmark
	    my ($fpr, $tpr) = $self->getBenchmarkRatesForOneScore($ratesAtEachScore, $score);

	    #add to $seqInfo hash for this peptide
	    $self->{Peptides}->{$type}->{$peptideSequenceId}->{$pOnePosition}->{"TPR"} = $tpr;
	    $self->{Peptides}->{$type}->{$peptideSequenceId}->{$pOnePosition}->{"FPR"} = $fpr;
	}
    }
}


############################################################################################################################################################
# getBenchmarkRatesForAllScores
# 
# Read benchmark score file that was output during a previous training run of the server to get TPRs and FPRs for each score.
# Benchmark score file has one line for each score (one score for each peptide that was trained on). Each line is of the format
# FP_Rate TP_Rate SVM_Score; these are tab-delimited)
#
# RETURN $scoreBins: hash containing results, of the format:
#                    $scoreBins->{$svmScore}->{tpRateAtScore} = $tpRate
#                    $scoreBins->{$svmScore}->{fpRateAtScore} = $fpRate 
#        $tpRate and $fpRate are both fractions
############################################################################################################################################################
sub getBenchmarkRatesForAllScores{
    my ($self) = @_;

    my $benchmarkScoreFile = $self->getParam("benchmark_score_file");
    my $benchmarkScoreFh = FileHandle->new("<" . $benchmarkScoreFile);

    unless ($benchmarkScoreFh =~ /FileHandle/){
	$self->writeError("file_missing", "getBenchmarkRatesForAllScores", "could not open benchmark score file name $benchmarkScoreFile: $!");
	return 0;
    }

    my $scoreBins;
    
    #read benchmark file that was output during training mode
    while (<$benchmarkScoreFh>){
	chomp;
	my $line = $_;
	next if ($line =~ /^\#/);

	my ($fpRate, $tpRate, $svmScore) = split('\t', $line);
	$svmScore = sprintf('%.3f', $svmScore);
	$scoreBins->{$svmScore}->{tpRate} = $tpRate;
	$scoreBins->{$svmScore}->{fpRate} = $fpRate;
    }
    return $scoreBins;
}


############################################################################################################################################################
# getBenchmarkRatesForOneScore
#
# Given an SVM score, return the FPR and TPR that it falls on in the benchmark set created by the server training mode.
# In a robustly trained server with lots of positives and negatives, there is a TPR and FPR for every score when rounded to the thousanth's place. Occasionally
# a score will be missing because no peptide in the training set received the score (this is especially the case at the high and low ends of the score distribution).
# In that case, take the TPR and FPR for the first score less than $svmScore. 
#
# If $svmScore is worse than anything in the benchmark set, just return the FPR and TPR for the worst score in the set (should be 1,1).
# 
# PARAM  $ratesAtEachScore: hash containing score bins (returned by $self->getBenchmarkRatesForAllScores())
# PARAM  $svmScore: score for which we are searching for benchmark values.
# RETURN (FPR, TPR) two element array.
############################################################################################################################################################
sub getBenchmarkRatesForOneScore{
    my ($self, $ratesAtEachScore, $svmScore) = @_;
 
    my @allScores = sort ({$a <=> $b} keys %$ratesAtEachScore);  #sorts from worst SVM score (most negative) to best
    my $worstScore = $allScores[0];
    #$svmScore = sprintf('%.3f', $svmScore);
    
    my $results = $ratesAtEachScore->{$svmScore};

    #if we don't have value for $svmScore, continue to decrement until we do
    while (!($results)){     
	$svmScore -= .001;
	$svmScore = sprintf('%.3f', $svmScore);   

	$results = $ratesAtEachScore->{$svmScore};
	if ($svmScore < $worstScore){ 
	    #we got something so bad that it is below the benchmark, just take the worst thing 
	    $results = $ratesAtEachScore->{$worstScore};
	}
    }

    my $tpRate = $results->{tpRate};
    my $fpRate = $results->{fpRate};
    
    return ($fpRate, $tpRate);
}


#######################################################################################################################################
# getPeptides()
# Returns all peptides with type $peptideType.
#
# PARAM  $peptideType: classification of peptides to return
# RETURN hash of the form $peptides->{$modbaseSeqId}->{$startPosition} = $seqInfo
#######################################################################################################################################

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


sub getGlobalErrorRegex{
    
    my ($self) = @_;
    return "file_missing|internal_error|invalid_model";
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

sub writeError{
    my ($self, $errorCode, $methodName, $errorMessage) = @_;
    $self->setErrors($errorCode);
    $self->writeLog("ERROR: SvmModel method $methodName Error: $errorMessage");
}


sub setErrors{
    my ($self, $errorKeyword) = @_;
    $self->{errors} = $errorKeyword;
}
sub getErrors{
    my ($self) = @_;
    my $errors = $self->{errors};
    return $errors;
}

#############################################################################################################################################################  
# Adds a sequence to the list of those for which no peptides were parsed in Application Scan mode (stored as hash).
# 
# PARAM  $modbaseSeqId: Sequence ID to add
# RETURN NULL
#############################################################################################################################################################  
sub writeNoPeptidesParsedForSeq{
    my ($self, $modbaseSeqId, $uniprotAccession) = @_;
    $self->{NoPeptidesParsedSeqs}->{$modbaseSeqId} = $uniprotAccession;
}

sub getNoPeptidesParsedSeqs{

    my ($self) = @_;
    my $noPeptidesParsedSeqs = $self->{NoPeptidesParsedSeqs};
    return $noPeptidesParsedSeqs;
}


sub trainsOnNegatives{
    return 1;
}

sub resetFeatureNumber{

    my ($self) = @_;

    $self->{MaxFeatureNumber} = 0;  
}


return 1;
