package PeptidePipeline;


use FastaReader;

use strict;
use FileHandle;
use POSIX qw(ceil floor);

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

    #error stuff
    $self->loadErrorCodes();

    #set global pipeline data structures
    $self->loadLog();
    $self->loadResultFh();

    $self->loadRunInfo();
}

#intended for application mode
sub parsePeptidesFromSequences{

    my ($self) = @_;

    $self->writeLog("Pipeline Application Mode; parsing peptides from protein sequences using rules file");
    my $startTime = time();

    my $peptideRulesFileName = $self->getParam("rules_file_name");
    my $peptideFastaFileName = $self->getParam("input_fasta_file_name");
    my $peptideLength = $self->getParam("peptide_length");
   

    my $peptideRules = $self->readPeptideRulesFile($peptideRulesFileName);

    $self->{PeptideRules} = $peptideRules;

    my $allPeptides = $self->parsePeptideFastaFile($peptideFastaFileName, $peptideLength, $peptideRules);

    $self->{AllPeptides} = $allPeptides;

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "parsePeptidesFromSequence");

}

#intended for training mode (update -- or application user-specified mode)
sub readPeptideInputFile{

    my ($self) = @_;

    $self->writeLog("Pipeline Training Mode; reading input file to retrieve user-specified peptides");
    my $startTime = time();
    
    #get input sequence info
    my $peptideInputFileName = $self->getParam("input_fasta_file_name");
    my $fastaReader = FastaReader->new($peptideInputFileName, 1);
    $self->handleReaderError($peptideInputFileName) unless $fastaReader;

    $fastaReader->read();
    my $sequences = $fastaReader->getSequences();

    my $keywordMismatch = $self->getParam("keyword_peptide_sequence_mismatch");

    my $allPeptides;

    #read sequences
    foreach my $sequenceInfo (keys %$sequences){
	my $sequenceList = $sequences->{$sequenceInfo};   #FastaReader either reads one line of full fasta sequence or multiple lines, each with part of sequence.  Stores results in array.
	                                                  #Here, it has read the whole line and stored it in one element array ($sequenceList)
	#parse sequence header
	my $sequence = $sequenceList->[0];
	my @cols = split('\|', $sequenceInfo);
	my $modbaseSeqId = $cols[0];
	
	my $uniprotAccession = $cols[1];
	my $peptides;
	my $peptideLength;

	#get all peptides to train on from this modbaseSeqId
	for (my $i = 2; $i < scalar(@cols); $i++){
	    my $nextPeptideSpec = $cols[$i];
	    my ($peptideStartPosition, $peptideSequence, $classification) = split("\_", $nextPeptideSpec);
	    next if ($peptideSequence eq $keywordMismatch);

	    $peptideLength = length($peptideSequence);
	    my $peptideEndPosition = $peptideStartPosition + $peptideLength - 1;

	    $peptides->{$peptideStartPosition}->{peptideSequence} = $peptideSequence;
	    $peptides->{$peptideStartPosition}->{peptideStartPosition} = $peptideStartPosition;   #this duplicates keys of $allPeptides but makes it easier to get the values when needed later
	    $peptides->{$peptideStartPosition}->{peptideEndPosition} = $peptideEndPosition;
	    $peptides->{$peptideStartPosition}->{classification} = $classification;
	}

	#create final $allPeptides data structure
	$allPeptides->{modbaseSeqId} = $modbaseSeqId;
	$allPeptides->{uniprotAccession} = $uniprotAccession;
	$allPeptides->{peptides} = $peptides;
	$allPeptides->{fullSequence} = $sequence;
	my $sequenceLength = length($sequence);
	$allPeptides->{sequenceLength} = $sequenceLength;
	$self->loadColumnInfo();
#	$self->createColumnInfo();  #putting this here for historical purposes.  TODO - see if it can go in initialize();
	$self->{AllPeptides} = $allPeptides;
    }

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "readPeptideInputFile");
}

sub printAllPeptides{
    my ($self) = @_;

    my $allPeptides = $self->{AllPeptides};

    my $modbaseSeq = $allPeptides->{modbaseSeqId};
    my $uniprotAccession = $allPeptides->{uniprotAccession};

    my $peptides = $allPeptides->{peptides};

    print STDERR "modbase: $modbaseSeq uniprot $uniprotAccession\n";

    foreach my $position (keys %$peptides){

	my $peptideInfo = $peptides->{$position};

	foreach my $keyName (keys %$peptideInfo){
	    my $value = $peptideInfo->{$keyName};
	    print STDERR "$keyName:\t$value\n";
	}
	print STDERR "END PEPTIDE\n";
    }
	
    

}

sub parsePeptideFastaFile{

    my ($self, $peptideFastaFile, $peptideLength, $peptideRules) = @_;
    
    my $fastaReader = FastaReader->new($peptideFastaFile, 1);
    $self->handleReaderError($peptideFastaFile) unless $fastaReader;
   
    $fastaReader->read();
    
    my $sequences = $fastaReader->getSequences();

    my $allPeptides;

    foreach my $sequenceInfo (keys %$sequences){
	my $sequenceList = $sequences->{$sequenceInfo};   #FastaReader either reads one line of full fasta sequence or multiple lines, each with part of sequence.  Stores results in array.
                                                          #Here, it has read the whole line and stored it in one element array ($sequenceList)
	my $sequence = $sequenceList->[0];
	my ($modbaseSeqId, $uniprotAccession) = split('\|', $sequenceInfo);
	
	my $peptides = $self->getPeptidesFromFullSequence($sequence, $peptideLength, $peptideRules);

	
	$allPeptides->{modbaseSeqId} = $modbaseSeqId;
	$allPeptides->{uniprotAccession} = $uniprotAccession;
	$allPeptides->{peptides} = $peptides;
	$allPeptides->{fullSequence} = $sequence;
	my $sequenceLength = length($sequence);
	$allPeptides->{sequenceLength} = $sequenceLength;
	
    }

    $self->loadColumnInfo();
#    $self->createColumnInfo();  #putting this here for historical purposes.  TODO - see if it can go in initialize();

    return $allPeptides;
}

sub getPeptidesFromFullSequence{
    
    my ($self, $sequence, $peptideLength, $peptideRules) = @_;
    
    my $peptides;

    my @seqArray = split('', $sequence);
    
    for (my $i = 0; $i < scalar(@seqArray) - $peptideLength; $i++){  #TODO -- validate that we will never have an index overrun with this OR just check the last peptide if there are no peptide rules and make sure we got it
	my $nextPeptide = substr($sequence, $i, $peptideLength);
	if ($self->passesCheck($nextPeptide, $peptideRules)){     #check to make sure peptide complies with rules
	    $peptides->{$i}->{peptideSequence} = $nextPeptide;
	    my $peptideEndPosition = $i + $peptideLength - 1;
	    $peptides->{$i}->{peptideStartPosition} = $i;   #this duplicates keys of $allPeptides but makes it easier to get the values when needed later
	    $peptides->{$i}->{peptideEndPosition} = $peptideEndPosition;
	    $peptides->{$i}->{classification} = "Application";
	}
    }
    return $peptides;
    
}

sub passesCheck{
    
    my ($self, $nextPeptide, $peptideRules) = @_;
    foreach my $positionNumber (keys %$peptideRules){
	my $bannedResidues = $peptideRules->{$positionNumber};
	my $residueAtPosition = substr($nextPeptide, $positionNumber - 1, 1);  #position is 1-based
	foreach my $bannedResidue (keys %$bannedResidues){
	    $bannedResidue = uc($bannedResidue);
	    if ($bannedResidue eq $residueAtPosition){
		return 0;
	    }
	}
    }
    return 1;
}

sub readPeptideRulesFile{
    my ($self,  $peptideRulesFileName) = @_;

    my $peptideRulesFh = FileHandle->new("<" . $peptideRulesFileName) || die "could not open peptide rules file $peptideRulesFileName\n"; #TODO -- better way to handle this and all other FileHandle die's 

    my $peptideRules;

    while (<$peptideRulesFh>){
	chomp;
	my $line = $_;               #line looks like 1 A C D F  where 1 is the position in the peptide to restrict and is followed by list of residues
	next if ($line =~ /^\#/);

	my @cols = split('\s+', $line);
	
	my $positionNumber = $cols[0];
	
	for (my $i = 1; $i < scalar(@cols); $i++){
	    my $nextResidue = $cols[$i];
	    $peptideRules->{$positionNumber}->{$nextResidue} = 1;
	}
    }
    return $peptideRules;

}


sub getBestModels{

    my ($self) = @_;
    
    $self->writeLog("Getting best models for sequences");
    my $startTime = time();

    my $modelTable = $self->getParam('model_table');

    my $fh = FileHandle->new("<" . $modelTable) || die "could not open model result table $modelTable\n";

    my $allPeptides = $self->{AllPeptides};
    my $peptides = $allPeptides->{peptides};
    my $modbaseSeqId = $allPeptides->{modbaseSeqId};
    my $sequenceLength = $allPeptides->{sequenceLength};
    my $peptideLength = $self->getParam("peptide_length");
    my $bestModelCriteria = $self->getParam("best_model_criteria");    #TODO - test other criteria
    
    my $modelInfo;
    
    #read through model file, get info for all models for this sequence
    my $totalModelCount = 0;
    while (<$fh>){
	chomp;
	my $line = $_;
	my @cols = split('\t', $line);
	my $seqId = $cols[0];
	if ($seqId eq $modbaseSeqId){
	    
	    $totalModelCount++;

	    my $alignmentId = $cols[1];
	    my $templateId = $cols[2];
	    my $modelId = $cols[3];
	    my $targetLength = $cols[4];
	    my $modelTargetStart = $cols[5];
	    my $modelTargetEnd = $cols[6];
	    my $sequenceIdentity = $cols[13];
	    my $eValue = $cols[14];
	    my $modelScore = $cols[15];
	    my $zDope = $cols[22];
	    my $nativeOverlap = $cols[25];
	    my $tsvmodMethod = $cols[26];

	    my $modelLength = $modelTargetEnd - $modelTargetStart;
	    my $coverage = ($modelLength * 1.0) / ($sequenceLength * 1.0);
	    my $roundedCoverage = sprintf("%.2f", $coverage);
	    
	    $modelInfo->{$modelId}->{alignmentId} = $alignmentId;
	    $modelInfo->{$modelId}->{templateId} = $templateId;
	    $modelInfo->{$modelId}->{modelTargetStart} = $modelTargetStart;
	    $modelInfo->{$modelId}->{modelTargetEnd} = $modelTargetEnd;
	    $modelInfo->{$modelId}->{templateSequenceIdentity} = $sequenceIdentity;
	    $modelInfo->{$modelId}->{eValue} = $eValue;
	    $modelInfo->{$modelId}->{modelScore} = $modelScore;
	    $modelInfo->{$modelId}->{zDope} = $zDope;
	    $modelInfo->{$modelId}->{nativeOverlap} = $nativeOverlap;
	    $modelInfo->{$modelId}->{tsvmodMethod} = $tsvmodMethod;

	    $modelInfo->{$modelId}->{modelLength} = $modelLength;
	    $modelInfo->{$modelId}->{coverage} = $roundedCoverage;
	}
    }
    $allPeptides->{totalModelCount} = $totalModelCount;
    
    #for each peptide, find best model that contains it
    foreach my $peptideStartPosition (keys %$peptides){
	my $peptideInfo = $peptides->{$peptideStartPosition};
	my $peptideEndPosition = $peptideInfo->{peptideEndPosition};
	my $modelsContainingPeptideCount;
       	my @modelsContainingPeptide;
	
	

	foreach my $modelId (keys %$modelInfo){
	    my $modelTargetStart = $modelInfo->{$modelId}->{modelTargetStart};
	    my $modelTargetEnd = $modelInfo->{$modelId}->{modelTargetEnd};
	
	    if ($modelTargetStart <= $peptideStartPosition && $modelTargetEnd >= $peptideEndPosition){  #we have found model containing peptide, add it to the list
		push(@modelsContainingPeptide, $modelId);
		$modelsContainingPeptideCount++;
	    }
	}

	#use model properties and critiera for best model (user specified) to get the model id of the highest quality model for this peptide
	my $bestModelId = $self->getBestModelByCriteria($bestModelCriteria, $modelInfo, @modelsContainingPeptide);
	my $modelTargetStart = $modelInfo->{$bestModelId}->{modelTargetStart};
	my $modelTargetEnd = $modelInfo->{$bestModelId}->{modelTargetEnd};
	#put all info for best model into $seqInfo
	if ($bestModelId){
	    foreach my $infoType (keys %{$modelInfo->{$bestModelId}}){
		$peptideInfo->{$infoType} = $modelInfo->{$bestModelId}->{$infoType};
	    }
	    $peptideInfo->{bestModelId} = $bestModelId;
	    
	    my $alignmentId = $peptideInfo->{alignmentId};
	    my $run = $self->getParam("modpipe_run_name");
	    
	    my $fullAlignmentFileName = $self->getAlignmentFileName($modbaseSeqId, $alignmentId, $run);
	    
	    my $templatePdbIds = $self->getTemplatePdbIds($modbaseSeqId, $peptideStartPosition, $fullAlignmentFileName, $run);
	    
	    $peptideInfo->{templatePdbIds} = $templatePdbIds;
	    $peptideInfo->{fullAlignmentFilePath} = $fullAlignmentFileName;
	    $peptideInfo->{modelsContainingPeptideCount} = $modelsContainingPeptideCount; 
	}

    }

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "getBestModels");
}


sub getAlignmentFileName{

    my ($self,  $modbaseSeqId, $alignmentId, $runName) = @_;

    my $runInfo = $self->getRunInfo();
    my $style = $runInfo->{$runName}->{style};
    my $path = $runInfo->{$runName}->{path};
       
    my $alignmentFileName;

    #retrieve full file path of alignment file; this will vary depending on which run style was used

    if ($style eq "old"){
	$alignmentId =~ /^(\S\S).*/;
	my $twoLetterCode = $1;

	$alignmentFileName = "$path" . "alignments/$twoLetterCode/$alignmentId/$alignmentId" . ".ali";
    }
    else {
	$modbaseSeqId =~ /^(\S\S\S).*/;
	my $threeLetterCode = $1;
	
	$alignmentFileName = "$path" . "$threeLetterCode/$modbaseSeqId/alignments/$alignmentId" . ".ali";
    }
    return $alignmentFileName;
}


sub getTemplatePdbIds{

    my ($self, $modbaseSeqId, $peptideStartPosition, $fullAlignmentFileName, $run) = @_;
    
    my $runInfo = $self->getRunInfo();
    my $runStyle = $runInfo->{$run}->{style};
    
    my ($targetInfoLine, $targetSequence, $templateInfoLine, $templateSequence) = $self->parseAlignmentFile($fullAlignmentFileName, $run);

    #new-style modpipe runs have target as first protein in ali file; old style has template; process file accordingly
    #figure out which line is the template info line
    my $lineToProcess;
    
    my @cols = split('\:', $templateInfoLine);
    my $pdbCode = $cols[1];
    my $proteinType = $cols[0];  #either 'structure' or 'structureX' depending on run type
    
    unless ($pdbCode =~ /^\S\S\S\S$/){    #expect 4 letter pdb code
	my $errorString = "Did not get four letter pdb code (got $pdbCode) in first column of pir line in expected position in alignment file $fullAlignmentFileName run style $runStyle";
	$self->writeError($modbaseSeqId, $peptideStartPosition, 0, "unexpected_pdb_code_in_alignment", "getBestModels", $errorString, 1);
	$pdbCode = "template_pdb_code_error";
    }
    unless ($proteinType =~ /structure/){   
	my $errorString = "Did not get protein type of 'structure' (got $proteinType) in first column of pir format in expected position in alignment file $fullAlignmentFileName run style $runStyle";
	$pdbCode = "template_pdb_code_error";
	$self->writeError($modbaseSeqId, $peptideStartPosition, 0, "incorrect_template_position_in_alignment", "getBestModels", $errorString, 1);
    }

    return $pdbCode;
}


sub getBestModelByCriteria{

    my ($self, $bestModelCriteria, $modelInfo, @modelsContainingPeptide) = @_;   #assumes the best is always the one with the greatest score; true for native overlap, coverage, model score, not DOPE (if we add it)

    my $currentBestCriteria = 0.0;
    my $currentBestModelId;
    foreach my $modelId (@modelsContainingPeptide){
	my $currentCriteria = $modelInfo->{$modelId}->{$bestModelCriteria};
	if ($currentCriteria > $currentBestCriteria){
	    $currentBestModelId = $modelId;
	    $currentBestCriteria = $currentCriteria;
	}
    }
    return $currentBestModelId;
}




sub parseDsspResults{
    my ($self) = @_;
    
    $self->writeLog("Parsing results of DSSP runs to get structure information from structures and models");
    my $startTime = time();

    my $dsspDir = $self->getParam("dssp_directory");   #relative to pipeline dir, ok to keep static for now

    my $peptides = $self->{AllPeptides}->{peptides};
    my $sequenceCount = 0;

    my $modbaseSeqId = $self->{AllPeptides}->{modbaseSeqId};    #TODO - decide if I want to parameterize this
    
    foreach my $peptideStartPosition (keys %$peptides){
	print STDERR "DSSP: getting position of next peptide $peptideStartPosition\n";

	my $peptideInfo = $peptides->{$peptideStartPosition};
	$sequenceCount++;
	$self->writeLog("Parsing DSSP results for sequence number $sequenceCount") if ($sequenceCount % 1000 == 0);
	#get dssp file and open
	$modbaseSeqId =~ /^(\S\S\S)\S+/;
	my $seqTopDir = $1;
	my $bestModelId = $peptideInfo->{bestModelId};
	
	next unless $bestModelId;
	next if $self->isError($modbaseSeqId, $peptideStartPosition, "no_model_in_expected_directory", $bestModelId);
	
	my $noFile = 0;
	my $dsspFile = "$dsspDir/$seqTopDir/$modbaseSeqId/$bestModelId.dssp";
	my $dsspFh = FileHandle->new("<" . $dsspFile) || ($noFile = 1);
	
	if ($noFile == 1){
	    print STDERR "could not open DSSP result file $dsspFile\n";   #this has happened occasionally and needs to be monitored.  Take out the die to do testing.
	    die;
	}
	#takes the annotated cleavage sequence start position and offsets if, if necessary, in accordance with modbase model residue numbering
	my ($peptideStartPosition, $peptideEndPosition) = $self->getModelCleavageRange($peptideInfo);
	$peptideStartPosition++;
	$peptideEndPosition++;     #convert from 0 to 1 based.  We might want to have it go in one based after parsing in step 1?
	my $peptideSequence = $peptideInfo->{peptideSequence};
	#move down past header info
	while (<$dsspFh>){
	    chomp;
	    my $line = $_;
	    last if ($line =~ /\s+\#/);
	}
	
	#read structure info
	my $residueSolventAcc = $self->getResidueSolventAcc();
	my $peptideSiteSSInfo;
	my $currentResiduePosition = 1; #starts at one no matter which style modpipe was used
	my $baseOnePeptideStart;  #peptide sites if pdb file were counted starting at one
	my $baseOnePeptideEnd;
	my $checkStructureString = "";
	my @structureValueArray;
	my @structureTypeArray;
	my @accessibilityFractionArray;
	my @accessibilityTypeArray;

	while (<$dsspFh>){
	    chomp;
	    my $line = $_;
	    
	    last if ($line =~ /^\s+$/);
	    
	    # 1135 1145   G  S >  S+     0   0   56  .*
	    
	    #count #res  #code   #type             #acc
	    if ($line =~ /\s+\d+\s+(\d+)\s+(\w)\s\s(.).{12}\s*\S+\s+(\d+).*/){
		
		my $residueNumber = $1;
		my $residueOneLetter = $2;
		my $structureType = $3;
		my $accessibility = $4;
		$currentResiduePosition++;

		$peptideSiteSSInfo->{$currentResiduePosition} = $structureType;
		
		#check if current residue is within the peptide sequence boundaries
		if ($residueNumber >= $peptideStartPosition && $residueNumber <= $peptideEndPosition){
		    
		    #quality control
		    my $positionInPeptideSeq = $residueNumber - $peptideStartPosition + 1;
		    $self->checkResidue($residueOneLetter, $positionInPeptideSeq, $peptideSequence, $modbaseSeqId, $peptideStartPosition);
		    $checkStructureString .= $structureType;
		    
		    #process secondary structure
		    my ($mappedStructureType, $mappedStructureValue) = $self->getMappedStructureType($structureType);
		    push (@structureTypeArray, $mappedStructureType);
		    push (@structureValueArray, $mappedStructureValue);

		    
		    #process solvent accessibility
		    my ($percentAccessible, $accCall) = $self->getSolventExposureFraction($accessibility, $residueOneLetter, $residueSolventAcc);
		    push (@accessibilityFractionArray, $percentAccessible);
		    push (@accessibilityTypeArray, $accCall);
		}
		
		#prepare for counting loop length later
		if ($residueNumber == $peptideStartPosition){
		    $baseOnePeptideStart = $currentResiduePosition;
		}
		if ($residueNumber == $peptideEndPosition){
		    $baseOnePeptideEnd = $currentResiduePosition;
		}
		
	    }
	    else {
		$self->writeLog("dssp file line is $line; regex did not work");
		die;
	    }
	}
	if ($checkStructureString eq ""){  #for now assume that if we couldn't get a result for secondary structure, the same is true for solvent accessibility
	    my $errorString = "Could not read DSSP secondary structure information from DSSP result file";
	    $self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId, "no_dssp_structure_info", "parseDsspResults", $errorString, 1);
	}
	
	#create strings from structure / accessibility values
	my $structureTypeString = join('', @structureTypeArray);
	my $structureValueString = join(',', @structureValueArray);
	my $accessibilityTypeString = join('', @accessibilityTypeArray);
	my $accessibilityFractionString = join(',', @accessibilityFractionArray);
	
	#set structure / accessibilty values for this cleavage sequence
	$peptideInfo->{structureTypes} = $structureTypeString;
	$peptideInfo->{structureValues} = $structureValueString;
	$peptideInfo->{accessibilityTypes} = $accessibilityTypeString;
	$peptideInfo->{accessibilityFractions} = $accessibilityFractionString;

#	my $loopCount = $self->countCleavageSiteLoopLength($peptideSiteSSInfo, $baseOnePeptideStart, $baseOnePeptideEnd);
	
	my $loopCount = "not counted";  #TODO - decide if we want to keep this
	$peptideInfo->{"loopLength"} = $loopCount;
    }

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "parseDsspResults");
}



sub getMappedStructureType{

    my ($self, $structureInputType) = @_;

    if ($structureInputType =~ /^\s+$/){
	return ("L", 1);
    }

    if ($structureInputType eq "H" || $structureInputType eq "G" || $structureInputType eq "I"){
	return ("A", 2);
    }
    if ($structureInputType eq "B" || $structureInputType eq "E"){
	return ("B", 3);
    }
    if ($structureInputType eq "S" || $structureInputType eq "T" || $structureInputType eq "L"){
	return ("L", 1);
    }
}

sub runPsipred{

    my ($self) = @_; 

    my $startTime = time();
    $self->writeLog("Running Psipred to get sequence-based secondary structure information");

    my $peptides = $self->{AllPeptides}->{peptides};
    my $sequenceCount = 0;

    my $modbaseSeqId = $self->{AllPeptides}->{modbaseSeqId};    #TODO - decide if I want to parameterize this

    my $psipredResultsDir = $self->getParam("psipred_results_dir");
    foreach my $peptideStartPosition (keys %$peptides){
	print STDERR "processing psipred peptide start position $peptideStartPosition\n";
	my $peptideInfo = $peptides->{$peptideStartPosition};
	
	my $sequenceCount = 0;
	
	my $peptideStartPosition = $peptideInfo->{peptideStartPosition};

	my $psipredStartPosition = $peptideStartPosition + 1;     #psipred expects 1-based residue numbering and ours go in as 0-based.   
	
	my $peptideSequence = $peptideInfo->{peptideSequence};
	my @peptideArray = split('', $peptideSequence);
	
	$sequenceCount++;
	$self->writeLog("Getting Psipred result for sequence number $sequenceCount") if ($sequenceCount % 10000 == 0);
	
	#Get psipred results file from modbase sequence id and input parameters
	$modbaseSeqId =~ /^(\S\S)\S+/;
	my $twoLetters = $1;
	my $fullPsipredDir = $psipredResultsDir . "/" . $twoLetters;
	my $fullPsipredFile = $fullPsipredDir . "/" . $modbaseSeqId . ".ss2";

	my $noFile = 0;
	
	if (-e $fullPsipredFile == 0){
	    $self->writeLog("Did not find Psipred file for modbase sequence ID $modbaseSeqId; running psipred command to create it");
	    $self->runPsipredCommand($modbaseSeqId, $fullPsipredDir);
	}
	
	my $psipredFh = FileHandle->new("<" . $fullPsipredFile) || ($noFile = 1);
	if ($noFile == 1){
	    $self->writeError($modbaseSeqId, $peptideStartPosition, 0, "no_psipred_result_file", "processPsipredResults", "Could not find Psipred Result File $fullPsipredFile", 1);
	    next;
	}


	my $peptideLength = $self->getParam("peptide_length");

	my $structureString = "";
	my @firstScoreArray;

	#scroll down to start of results
	for (my $i = 0; $i < 2; $i++){
	    <$psipredFh>;
	}
	
	while (<$psipredFh>){
	    chomp;
	    my $line = $_;
	    if ($line =~ /^\s+(\S.*)/){
		$line = $1;   #trim initial blank space form line for three digit residue numbers (four digit not an issue)
	    }
	    
	    my ($residueNumber, $residueOneLetter, $structureCall, $firstScore, $secondScore, $thirdScore)  = split('\s+', $line);  
			    
	    if ($residueNumber == $psipredStartPosition){  #found cleavage sequence position, process the entire thing within the next for loop
		for (my $i = 0; $i < $peptideLength; $i++){
		    my $expectedPeptideResidue = $peptideArray[$i];
		    if ($expectedPeptideResidue ne $residueOneLetter){
			my $errorString =  "Error: psipred cleavage sequence letter $residueOneLetter position $residueNumber does not equal expected cleavage seq entry\n";
			$errorString .= "$expectedPeptideResidue in $peptideSequence position $i in modbase sequence $modbaseSeqId\n";
			$structureString = "psipred sequence error";  #fix this according to notes in comments at bottom of file
			
			$self->writeError($modbaseSeqId, $peptideStartPosition, 0, "psipred_sequence_mismatch", "processPsipredResults", $errorString, 1);
			@firstScoreArray = ();
			last;
			
		    }
		    elsif ($structureCall eq "C"){
			$structureString .= "L";
		    }
		    elsif ($structureCall eq "E"){
			$structureString .= "B";
		    }
		    elsif ($structureCall eq "H"){
			$structureString .= "A";
		    }
		    else {
			
			die "did not get expected psipred call for modbase seq id $modbaseSeqId position $residueNumber (should either be 'C', 'E', or 'H')\n";
		    }
		    
		    push(@firstScoreArray, $firstScore);
		    $line = <$psipredFh>;
		    chomp $line;
		    if ($line =~ /^\s+(\S.*)/){
			$line = $1;   #trim initial blank space form line for three digit residue numbers (four digit not an issue)
		    }
		    ($residueNumber, $residueOneLetter, $structureCall, $firstScore, $secondScore)  = split('\s+', $line); #read in next line
		}
		last; #we have finished the cleavage sequence.  Theoretically if the protein has more than one cleavage site we could keep going but will just process this sequence again instead
	    }
	}
	if ($structureString eq ""){   #did not find expected cleavage seq start position in psipred result filexo
	    $self->writeError($modbaseSeqId, $peptideStartPosition, 0, "no_psipred_cleavage_sequence", "processPsipredResults", "Found Psipred results file but did not find cleavage sequence start position $peptideStartPosition in file", 1); 
	    
	}
	else {
	    my $firstScoreString = join (',', @firstScoreArray);
	    
	    $peptideInfo->{psipredString} = $structureString;
	    
	    $peptideInfo->{psipredFirstScores} = $firstScoreString;
	    #will probably be adding psipred score soon as well
	}
    }
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "runPsipred");

}


sub runDisopred {

    my ($self) = @_; 

    my $startTime = time();
    $self->writeLog("Running disopred to get sequence-based disorder information");

    my $peptides = $self->{AllPeptides}->{peptides};
    my $sequenceCount = 0;

    my $modbaseSeqId = $self->{AllPeptides}->{modbaseSeqId};    #TODO - decide if I want to parameterize this

    my $disopredResultsDir = $self->getParam("disopred_results_dir");
    foreach my $peptideStartPosition (keys %$peptides){
	my $peptideInfo = $peptides->{$peptideStartPosition};

	my $sequenceCount = 0;
	
	my $peptideStartPosition = $peptideInfo->{peptideStartPosition};
	my $disopredStartPosition = $peptideStartPosition + 1;     #disopred expects 1-based residue numbering and ours go in as 0-based.   TODO -- see if we want to change this.

	my $peptideSequence = $peptideInfo->{peptideSequence};
	my @peptideArray = split('', $peptideSequence);
	
	$sequenceCount++;
	$self->writeLog("Getting Disopred result for sequence number $sequenceCount") if ($sequenceCount % 10000 == 0);

	#Get disopred results file from modbase sequence id and input parameters
	$modbaseSeqId =~ /^(\S\S)\S+/;
	my $twoLetters = $1;
	my $fullDisopredDir = $disopredResultsDir . "/" . $twoLetters;
	my $fullDisopredFile = $fullDisopredDir . "/" . $modbaseSeqId . ".diso";

	my $noFile = 0;
	
	if (-e $fullDisopredFile == 0){
	    $self->writeLog("Did not find Disopred file for modbase sequence ID $modbaseSeqId; running disopred command to create it");
	    $self->runDisopredCommand($modbaseSeqId, $fullDisopredDir);
	}

	my $disopredFh = FileHandle->new("<" . $fullDisopredFile) || ($noFile = 1);
	if ($noFile == 1){
	    $self->writeError($modbaseSeqId, $peptideStartPosition, 0, "no_disopred_result_file", "processDisopredResults", "Could not find Disopred Result File $fullDisopredFile", 1);
	    next;
	}

	my $peptideLength = $self->getParam("peptide_length");

	my $disorderString = "";
	my @firstScoreList;
	#scroll down to start of results
	for (my $i = 0; $i < 5; $i++){
	    <$disopredFh>;
	}
	

	while (<$disopredFh>){
	    chomp;
	    my $line = $_;
	    my ($blank, $residueNumber, $residueOneLetter, $disorderCall, $firstScore, $secondScore)  = split('\s+', $line);  #one-based residueNumber
	    
	    if ($residueNumber == $disopredStartPosition){  #found cleavage sequence position, process the entire thing within the next for loop
		
		for (my $i = 0; $i < $peptideLength; $i++){
		    my $expectedPeptideResidue = $peptideArray[$i];
		    if ($expectedPeptideResidue ne $residueOneLetter){

			my $errorString =  "Error: disopred peptide letter $residueOneLetter position $residueNumber does not equal expected peptide entry\n";
			$errorString .= "$expectedPeptideResidue in $peptideSequence position $i in modbase sequence $modbaseSeqId\n";
			$disorderString = "disopred sequence error";  #fix this according to notes in comments at bottom of file
			$self->writeError($modbaseSeqId, $peptideStartPosition, 0, "disopred_sequence_mismatch", "processDisopredResults", $errorString, 1);
			@firstScoreList = ();
			last;
		    }

		    elsif ($disorderCall eq "*"){
			$disorderString .= "D";
		    }
		    elsif ($disorderCall eq "."){
			$disorderString .= "O";
			}
		    else {
			
			die "did not get expected disopred call for modbase seq id $modbaseSeqId position $residueNumber (should either be '*' or '.')\n";  #TODO - make into error code
		    }
		    push (@firstScoreList, $firstScore);
		    $line = <$disopredFh>;
		    chomp $line;
		    ($blank, $residueNumber, $residueOneLetter, $disorderCall, $firstScore, $secondScore)  = split('\s+', $line); #read in next line
		}
		last; #we have finished the cleavage sequence.  Theoretically if the protein has more than one cleavage site we could keep going but will just process this sequence again instead
	    }
	}
	if ($disorderString eq ""){   #did not find expected cleavage seq start position in disopred result file
	
	    $peptideInfo->{disopredString} = "disopred no result error";
	}
	else {
	    $peptideInfo->{disopredString} = $disorderString;

	    my $firstScoreString = join(',', @firstScoreList);
	    $peptideInfo->{disopredFirstScores} = $firstScoreString;
	    #will probably be adding disopred score soon as well
	}
	    
    }
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "runDisopred");
}

sub runDisopredCommand{
    my ($self, $seqId, $outputDir) = @_;
    

    my $runDisopredCmd = $self->getParam("run_disopred_cmd");
    my $peptideFastaFileName = $self->getParam("input_fasta_file_name");


    #disopred output files have the same prefix as the fasta file name prefix.  TODO - check this works if we have a fasta file name without .fasta extension
    my $shortSeqId = $peptideFastaFileName;
    if ($peptideFastaFileName =~ /(.*)\.fasta/){
	$shortSeqId = $1;
    }

    my $cmd = "$runDisopredCmd $peptideFastaFileName";

    $self->writeLog("running disopred command $cmd");

    system($cmd);

    my $mkdirCmd = "mkdir -p $outputDir";
    $self->writeLog("making new disopred output directory $mkdirCmd");
    system($mkdirCmd);
    
    my $horizCpCmd = "cp $shortSeqId.horiz_d $outputDir/$seqId.horiz_d";
    $self->writeLog("copy disopred horiz_d results with $horizCpCmd");
    system($horizCpCmd);

    my $disoCpCmd = "cp $shortSeqId.diso $outputDir/$seqId.diso";
    $self->writeLog("copy disopred diso results with $disoCpCmd");
    system($disoCpCmd);

}



sub runPsipredCommand{
    my ($self, $seqId, $outputDir) = @_;

    my $runPsipredCmd = $self->getParam("run_psipred_cmd");
    my $peptideFastaFileName = $self->getParam("input_fasta_file_name");


    #psipred output files have the same prefix as the fasta file name prefix.  TODO - check this works if we have a fasta file name without .fasta extension
    my $shortSeqId = $peptideFastaFileName;
    if ($peptideFastaFileName =~ /(.*)\.fasta/){
	$shortSeqId = $1;
    }

    my $cmd = "$runPsipredCmd $peptideFastaFileName";

    $self->writeLog("running psipred command $cmd");

    system($cmd);

    my $mkdirCmd = "mkdir -p $outputDir";
    $self->writeLog("making new psipred output directory $mkdirCmd");
    system($mkdirCmd);

    my $horizCpCmd = "cp $shortSeqId.horiz $outputDir/$seqId.horiz";
    $self->writeLog("copy psipred horiz results with $horizCpCmd");
    system($horizCpCmd);
   

    my $ss2CpCmd = "cp $shortSeqId.ss2 $outputDir/$seqId.ss2";
    $self->writeLog("copy psipred ss2 results with $ss2CpCmd");
    system($ss2CpCmd);

}



sub applyGoAnnotation{


}


sub makeColumnHeaderString{

    my ($self) = @_;

    my $columnHeaderString = "";
    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");

    my $columns = $self->{ColumnInfo};
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   #check if the output file for Peptide Pipeline mode is supposed to output this column type
	if ($modes =~ /$keywordFeatureInternal/){
	    my $displayName = $columns->{$columnShortName}->{displayName};
	    my $method = $columns->{$columnShortName}->{method};
	    if ($method eq "universal"){
		$columnHeaderString .= $displayName . "\t";
	    }
	    else {
		if ($self->runMethod($method, 1)){
		    $columnHeaderString .= $displayName . "\t";
		}
	    }
	}
    }
    return $columnHeaderString;
}


sub makeOutputLine{
    my ($self, $outputInfo) = @_;
    my $columns = $self->{ColumnInfo};
    my $outputLine = "";
    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   #check if the output file for Peptide Pipeline mode is supposed to output this column type
	if ($modes =~ /$keywordFeatureInternal/){
	    my $method = $columns->{$columnShortName}->{method};
	
	    if ($method eq "universal"){
	    
		my $nextOutputValue = $outputInfo->{$columnShortName};
		$outputLine .= $nextOutputValue . "\t";
	    }
	    else {
		if ($self->runMethod($method, 1)){
		    my $nextOutputValue = $outputInfo->{$columnShortName};
		    $outputLine .= $nextOutputValue . "\t";
		}
	    }
	}
    }
    return $outputLine;
}


#taken from Rose et al, Science, 1985, which determined these values using atomic radii calculated in Lee & Richards, 1971
sub getResidueSolventAcc{

    my ($self) = @_;

    my $residues;

    $residues->{"A"} = 118.1;
    $residues->{"C"} = 146.1;
    $residues->{"D"} = 158.7;
    $residues->{"E"} = 186.2;

    $residues->{"F"} = 222.8;
    $residues->{"G"} = 88.1;
    $residues->{"H"} = 202.5;
    $residues->{"I"} = 181.0;

    $residues->{"K"} = 225.8;
    $residues->{"L"} = 193.1;
    $residues->{"M"} = 203.4;
    $residues->{"N"} = 165.5;

    $residues->{"P"} = 146.8;
    $residues->{"Q"} = 193.2;
    $residues->{"R"} = 256.0;
    $residues->{"S"} = 129.8;

    $residues->{"T"} = 152.5;
    $residues->{"V"} = 164.5;
    $residues->{"W"} = 266.3;
    $residues->{"Y"} = 236.8;

    return $residues;

}

sub getSolventExposureFraction{
    my ($self, $accessibility, $residueName, $residueSolventAcc) = @_;

    $accessibility *= 1.0;
		
    my $totalSa = $residueSolventAcc->{$residueName};
    my $percentSa = $accessibility / $totalSa;
       
    my $accCall;
    if ($percentSa > 0.33){
	$accCall .= "A";
    }
    else {
	$accCall .= "N";
    }

    return ($percentSa, $accCall);
}




sub writeUserResults{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Writing results");

    return unless $self->runMethod("writeUserResults");
    
    my $resultsFh = $self->{ResultFh};
    
    my $columnHeaderString = $self->makeColumnHeaderString();
    print $resultsFh $columnHeaderString . "\n";

    my $cleavageSeqLength = $self->{CleavageSequenceLength};
    my $residueSolventAcc = $self->getResidueSolventAcc();

    my $peptides = $self->{AllPeptides}->{peptides};
    my $sequenceCount = 0;

    my $modbaseSeqId = $self->{AllPeptides}->{modbaseSeqId};    #TODO - decide if I want to parameterize this
    my $uniprotAccession = $self->{AllPeptides}->{uniprotAccession};
    my $totalModelCount = $self->{AllPeptides}->{totalModelCount};

    my $disopredResultsDir = $self->getParam("disopred_results_dir");  #TODO -- check if we are running disopred before asking for this param
    foreach my $peptideStartPosition (keys %$peptides){
	
	my $peptideInfo = $peptides->{$peptideStartPosition};
	
	my $bestModelId = $peptideInfo->{bestModelId};
	
	my $outputInfo;
	
	#general sequence information    
	my @attributes = ("alternateAccession", "spliceAccession", "modelsContainingPeptideCount", "run",  "modelTargetStart", "modelTargetEnd", "peptideStartPosition", "peptideEndPosition", "peptideSequence",  "modelScore", "templateSequenceIdentity", "coverage", "loopLength", "sharedProteaseInfo", "classification",  "description",  "disopredString", "disopredFirstScores",  "psipredString", "psipredFirstScores",  "fullAlignmentFilePath");
	    
	foreach my $attribute (@attributes){
	    my $value = $peptideInfo->{$attribute};
	    $outputInfo = $self->addOutputValue($outputInfo, $value, $attribute);
	}
	my $url = "";
	if ($bestModelId){
	    $url = "http://salilab.org/modbase/search?modelID=" . $bestModelId . "&displaymode=moddetail";
	}
	
	$outputInfo = $self->addOutputValue($outputInfo, $modbaseSeqId, "sequenceId");
	$outputInfo = $self->addOutputValue($outputInfo, $uniprotAccession, "uniprotAccession");
	$outputInfo = $self->addOutputValue($outputInfo, $totalModelCount, "totalModelCount");
	$outputInfo = $self->addOutputValue($outputInfo, $bestModelId, "modelId");
	$outputInfo = $self->addOutputValue($outputInfo, $url, "url");
	
	
	if ($bestModelId){
	    
	    my @structureAttributes = ("zDope", "nativeOverlap", "tsvmodMethod", "templatePdbIds", "similarityString", "templateCleavageSequence", "structureTypes", "accessibilityTypes", "structureValues", "accessibilityFractions");
		
	    foreach my $structureAttribute (@structureAttributes){
		my $value = $peptideInfo->{$structureAttribute};
		$outputInfo = $self->addOutputValue($outputInfo, $value, $structureAttribute);
	    }
	}
	
	#cell lines
	if ($self->getParam("addLysateExpression") eq "yes"){
	    my $cellLinesUsed = $self->getParam("cell_lines_to_use");
	    my @cellLineNames = split(',', $cellLinesUsed);
	    
	    my $cellLineInfo = $peptideInfo->{cellLines};
	    foreach my $cellLineName (@cellLineNames){   #the order will be correct because both the column names and the values were determined by the cell_lines_to_use param
		my $callString = "";
		my $affyIdString = "";
		my $pCount = 0;
		my $affyIds = $cellLineInfo->{$cellLineName};
		foreach my $affyId(sort ({$affyIds->{$b} cmp $affyIds->{$a}} keys %$affyIds)){
		    my $call = $affyIds->{$affyId};
		    $callString .= $call . " ";
		    $affyIdString .= $affyId . " ";
		    if ($call eq 'P'){
			$pCount++;
		    }
		}
		my $finalString = "($pCount) . " . $callString . ": " . $affyIdString;
		$outputInfo = $self->addOutputValue($outputInfo, $finalString, $cellLineName); 
	    }
	}
	
	#flags
	my $flagNames = $self->{FlagNames};
	my $outputFlags = $self->{OutputFlags};
	foreach my $flagNumber (sort keys %$flagNames){
	    my $flagName = $flagNames->{$flagNumber};
	    if ($outputFlags->{$flagName} eq "yes"){
		
		my $flagValue = $peptideInfo->{$flagName};
		$outputInfo = $self->addOutputValue($outputInfo, $flagValue, $flagName);
	    }
	}

	my $errors = $peptideInfo->{errors};
	my $errorString = "";
	if ($errors){
	    $errorString = join (',', @$errors);
	}
	else {
	    
	    $errorString = $self->getParam("keyword_no_cluster_errors");
	}

	$outputInfo = $self->addOutputValue($outputInfo, $errorString, "errors");
	my $outputLine = $self->makeOutputLine($outputInfo);
	print $resultsFh $outputLine . "\n";


    }
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "writeUserResults");
}



sub createColumnInfo{
    
    my ($self) = @_;

    my $columnInfo;
    
    $self->addColumn("uniprotAccession", "Uniprot Accession", "normal");
    $self->addColumn("alternateAccession", "Alternate ID", "normal");
    $self->addColumn("spliceAccession", "Sequence Splice Form", "normal");
    $self->addColumn("totalModelCount", "Total Models For Sequence", "normal");
    $self->addColumn("modelsContainingPeptideCount", "Models Containing Peptide", "normal");
    $self->addColumn("run", "Dataset", "normal");
    $self->addColumn("modelTargetStart", "Target Start", "normal");
    $self->addColumn("modelTargetEnd", "Target End", "normal");
    $self->addColumn("peptideStartPosition", "Peptide Start", "normal");
    $self->addColumn("peptideEndPosition", "Peptide End", "normal");
    $self->addColumn("peptideSequence", "Peptide Sequence", "normal");
    $self->addColumn("hmmScore", "HMM Score", "applyHmm");
    $self->addColumn("hmmEValue", "HMM E-Value", "applyHmm");
    $self->addColumn("modelScore", "Model Score", "normal");
    $self->addColumn("zDope", "Dope Score", "normal");
    $self->addColumn("nativeOverlap", "TSVMod Native Overlap", "normal");
    $self->addColumn("tsvmodMethod", "TSVMod Method", "normal");
    $self->addColumn("templateSequenceIdentity", "Template Sequence Identity", "normal");
    $self->addColumn("coverage", "Model Coverage", "normal");
    $self->addColumn("loopLength", "Loop Length", "normal");
    $self->addColumn("sharedProteaseInfo", "Shares Site With Other Protease", "normal");
    $self->addColumn("description", "Name", "normal");
    $self->addColumn("templatePdbIds", "Template PDB ID", "normal");
    $self->addColumn("similarityString", "Peptide Similarity To Template", "normal");
    $self->addColumn("templateCleavageSequence", "Corresponding Sequence in Template", "normal");
    $self->addColumn("structureTypes", "Peptide Secondary Structure Type", "runDssp");
    $self->addColumn("accessibilityTypes", "Peptide Accessibility Prediction", "runDssp");
    $self->addColumn("structureValues", "Peptide Structure Values", "runDssp");
    $self->addColumn("accessibilityFractions", "Peptide Predicted Accessibility Fraction", "runDssp");
    

    if ($self->getParam("addLysateExpression") eq "yes"){
	my $cellLinesUsed = $self->getParam("cell_lines_to_use");
	my @cellLineColumnNames = split(',', $cellLinesUsed);
	foreach my $cellLine (@cellLineColumnNames){
	    $self->addColumn($cellLine, $cellLine, "addLysateExpression");
	}
    }

    $self->addColumn("disopredString", "Disopred Prediction", "processDisopredResults");
    $self->addColumn("disopredFirstScores", "Disopred Scores", "processDisopredResults");
    $self->addColumn("profAccString", "Prof Three-State Accessibility", "processProfResults");
    $self->addColumn("relativeAccScore", "Prof Scores", "processProfResults");
    $self->addColumn("profReliability", "Prof Reliability", "processProfResults");
    $self->addColumn("psipredString", "PSIPRED Prediction", "processPsipredResults");
    $self->addColumn("psipredFirstScores", "PSIPRED Scores", "processPsipredResults");
    $self->addColumn("availableAntibodyName", "Antibody Name", "addAntibodyAvailability");
    $self->addColumn("availableAntibodyProductNumber", "Antibody ProductNumber", "addAntibodyAvailability");
    $self->addColumn("domainBoundaries", "Domain Boundaries","processDomainAnnotation");
    $self->addColumn("domainBoundaryValue", "Boundary Cleavage Value", "processDomainAnnotation");
    $self->addColumn("goIdString", "GO IDs", "processGoAnnotation");
    $self->addColumn("goDescriptionString", "GO Functions", "processGoAnnotation");
    $self->addColumn("classification", "Classification", "normal");
    
    my $flagNames = $self->{FlagNames};
    my $outputFlags = $self->{OutputFlags};

    foreach my $flagNumber (sort keys %$flagNames){
	my $flagName = $flagNames->{$flagNumber};
	if ($outputFlags->{$flagName} eq "yes"){
	    $self->addColumn($flagName, $flagName, "flag");
	}
    }

    $self->addColumn("fullAlignmentFilePath", "Alignment File Path", "normal");
    $self->addColumn("sequenceId", "Sequence ID", "normal");
    $self->addColumn("modelId", "Model ID", "normal");
    $self->addColumn("url", "Model URL", "normal");

}

sub addColumn{

    my ($self, $internalId, $display, $type) = @_;

    if (lc($display) eq "status"){
	die "cannot add the column 'status' to output; it is a reserved word for benchmarking\n";   #until we make classes for all of these attributes this is the best I can do
    }

    #check to make sure column output name is valid from column info
    #set display order to be read from line

    $self->{Columns}->{$internalId}->{display} = $display;
    $self->{Columns}->{$internalId}->{type} = $type;

    my $currentCounter = $self->{ColumnCounter};
    $self->{Columns}->{$internalId}->{displayOrder} = $currentCounter;
    $self->{ColumnCounter}++;
}

sub addOutputValue{
    my ($self, $outputInfo, $value, $internalName) = @_;
    
    my $columns = $self->{ColumnInfo};
    my $columnInfo = $columns->{$internalName};
    unless ($columnInfo){
	my $keyString = join(", ", keys %$columns);
	print STDERR "ERROR:  trying to add value to column info with internal id $internalName which is not a defined internal id.  List of these IDs:\n$keyString\n";
	exit(1);
    }
    $outputInfo->{$internalName} = $value;

    return $outputInfo;
}


sub parseAlignmentFile{

    my ($self, $fullAlignmentFileName, $runName) = @_;

    
    my $runInfo = $self->getRunInfo();
    my $runStyle = $runInfo->{$runName}->{style};
    
    my $alignmentFh = FileHandle->new("<" . $fullAlignmentFileName) || die "could not open $fullAlignmentFileName\n";
    
    
    #read in first protein
    my $blank = <$alignmentFh>;
    my $firstHeaderLine = <$alignmentFh>;
    my $firstInfoLine = <$alignmentFh>;
    my $firstSequence = "";

    #read in protein sequence
    while (<$alignmentFh>){ 
	my $line = $_;
	chomp $line;
	$firstSequence .= $line;
	if ($line =~ /\*/){
	    $firstSequence =~ /(.*)\*$/;  #trim trailing * from sequence we have already built
	    $firstSequence = $1;
	    last;
	}	
    }

    #read in second protein info
    my $blank = <$alignmentFh>;
    my $secondHeaderLine = <$alignmentFh>;
    my $secondInfoLine = <$alignmentFh>;
    my $secondSequence = "";

    #read in protein sequence
    while (<$alignmentFh>){
	my $line = $_;
	chomp $line;
	$secondSequence .= $line;
	if ($line =~ /\*/){
	    $secondSequence =~ /(.*)\*$/;  #trim trailing *
	    $secondSequence = $1;
	    last;
	}
    }

    my ($templateInfoLine, $targetInfoLine, $templateSequence, $targetSequence);
    
    if ($runStyle eq "old"){
	$templateInfoLine =  $firstInfoLine;
	$templateSequence = $firstSequence;

	$targetInfoLine =  $secondInfoLine;
	$targetSequence = $secondSequence;
    }
    elsif ($runStyle eq "new"){
	$targetInfoLine =  $firstInfoLine;
	$targetSequence = $firstSequence;

	$templateInfoLine =  $secondInfoLine;
	$templateSequence = $secondSequence;
    }
    else {
	die "did not get expected run style 'old' or 'new' for run $runName\n";
    }


    return ($targetInfoLine, $targetSequence, $templateInfoLine, $templateSequence);


}



#taken from Rose et al, Science, 1985, which determined these values using atomic radii calculated in Lee & Richards, 1971
sub getResidueSolventAcc{

    my ($self) = @_;

    my $residues;

    $residues->{"A"} = 118.1;
    $residues->{"C"} = 146.1;
    $residues->{"D"} = 158.7;
    $residues->{"E"} = 186.2;

    $residues->{"F"} = 222.8;
    $residues->{"G"} = 88.1;
    $residues->{"H"} = 202.5;
    $residues->{"I"} = 181.0;

    $residues->{"K"} = 225.8;
    $residues->{"L"} = 193.1;
    $residues->{"M"} = 203.4;
    $residues->{"N"} = 165.5;

    $residues->{"P"} = 146.8;
    $residues->{"Q"} = 193.2;
    $residues->{"R"} = 256.0;
    $residues->{"S"} = 129.8;

    $residues->{"T"} = 152.5;
    $residues->{"V"} = 164.5;
    $residues->{"W"} = 266.3;
    $residues->{"Y"} = 236.8;

    return $residues;

}

sub getSolventExposureFraction{
    my ($self, $accessibility, $residueName, $residueSolventAcc) = @_;

    $accessibility *= 1.0;
		
    my $totalSa = $residueSolventAcc->{$residueName};
    my $percentSa = $accessibility / $totalSa;
       
    my $accCall;
    if ($percentSa > 0.33){
	$accCall .= "A";
    }
    else {
	$accCall .= "N";
    }

    return ($percentSa, $accCall);
}

sub getMappedStructureType{

    my ($self, $structureInputType) = @_;

    if ($structureInputType =~ /^\s+$/){
	return ("L", 1);
    }

    if ($structureInputType eq "H" || $structureInputType eq "G" || $structureInputType eq "I"){
	return ("A", 2);
    }
    if ($structureInputType eq "B" || $structureInputType eq "E"){
	return ("B", 3);
    }
    if ($structureInputType eq "S" || $structureInputType eq "T" || $structureInputType eq "L"){
	return ("L", 1);
    }
}



sub getModelOffset{
 
    my ($self, $peptideInfo) = @_;

    my $run = $self->getParam("modpipe_run_name");    #TODO -- currently this is parameterized, but might need to be dynamic if we add more runs and get best models accordingly
    my $runInfo = $self->getRunInfo();
    if (($run ne "human") && ($run ne "human_4-2007") && ($run ne "snp-human2")  && ($run ne "snp-human3")  && ($run ne "snp-human4") && ($run ne "trembl2004") && ($run ne "human_2008" )){
	die "run for this sequence not of expected type: $run\n";
    }
    
    my $hasOffset = $runInfo->{$run}->{offset};
    
    my $offset;

    if ($hasOffset){
	$offset = $peptideInfo->{modelTargetStart} - 1;
    }
    else {
	$offset = 0;
    }
    return $offset;
}




sub getModelCleavageRange{

    my ($self, $peptideInfo) = @_;

    my $offsets;

    my $offset = $self->getModelOffset($peptideInfo);
    
    my $peptideStartPosition = $peptideInfo->{peptideStartPosition};
    my $peptideEndPosition = $peptideInfo->{peptideEndPosition};
    
    $peptideStartPosition -= $offset;
    $peptideEndPosition -= $offset;

    return ($peptideStartPosition, $peptideEndPosition);
}

sub checkInput{
    my ($self) = @_;
    my $allPeptides = $self->{AllPeptides};

    my $peptides = $allPeptides->{peptides};
    my $peptideCount = scalar(keys %$peptides);
    if ($peptideCount < 1){
	return 0;
    }
    else {
	return 1;
    }
}

sub writeNoInput{
    my ($self) = @_;

    my $resultsFh = $self->{ResultFh};
    
    my $columnHeaderString = $self->makeColumnHeaderString();
    print $resultsFh $columnHeaderString . "\n";

    my $noPeptidesParsedKeyword = $self->getParam("keyword_no_peptides_parsed");
    print $resultsFh $noPeptidesParsedKeyword . "\n";
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
    my $logFileName =$self->getParam("peptide_pipeline_log_name");
    my $fullLogFileName = $pipelineDir . "/" . $run . "/" . $logFileName;
 
    my $logFh = FileHandle->new(">" . $fullLogFileName) || die "could not open $fullLogFileName for writing\n";
    $self->{LogFileName} = $fullLogFileName;
    $self->{Log} = $logFh;
}

sub loadResultFh{
    my ($self) = @_;

    my $run = $self->getParam("run_name");
    my $pipelineDir = $self->getParam("pipeline_directory");
    my $resultsFileName = $self->getParam("peptide_pipeline_result_file_name");
    my $fullResultsFileName = $pipelineDir . "/" . $run . "/" . $resultsFileName;
    my $resultsFh = FileHandle->new(">" . $fullResultsFileName) || die "could not open $fullResultsFileName for writing\n";
    $self->{ResultsFileName} = $fullResultsFileName;
    $self->{ResultFh} = $resultsFh;

}



sub loadErrorCodes{

    my ($self) = @_;
    $self->{ErrorCodes}->{"no_dssp_structure_info"} = "Sequences for which there were models but DSSP secondary structure output file could not be read at the cleavage sequence position.";
    $self->{ErrorCodes}->{"no_model_in_expected_directory"} = "Sequences for which there was a best model but it was not found in expected source modbase directory.";

    $self->{ErrorCodes}->{"no_disopred_result_file"} = "Sequences for which no Disopred results file was found.";
    $self->{ErrorCodes}->{"disopred_sequence_mismatch"} = "Sequences for which the residues at the cleavage sequence positions in the Disopred results file did not match up with the annotated cleavage sequence.";

    $self->{ErrorCodes}->{"no_prof_result_file"} = "Sequences for which no Prof results file was found.";
    $self->{ErrorCodes}->{"prof_sequence_mismatch"} = "Sequences for which the residues at the cleavage sequence positions in the Prof results file did not match up with the annotated cleavage sequence.";

    $self->{ErrorCodes}->{"no_psipred_result_file"} = "Sequences for which no Psipred results file was found.";
    $self->{ErrorCodes}->{"psipred_sequence_mismatch"} = "Sequences for which the residues at the cleavage sequence positions in the Psipred results file did not match up with the annotated cleavage sequence.";
    $self->{ErrorCodes}->{"no_psipred_cleavage_sequence"} = "Sequences for which there was a PSIPRED results file but the start position of the cleavage sequence was not found in the file.";

    $self->{ErrorCodes}->{"unexpected_pdb_code_in_alignment"} = "Alignment files that did not have a four letter pdb code in the expected column in the pir line for the template.";
    $self->{ErrorCodes}->{"incorrect_template_position_in_alignment"} = "Alignment files that appear to have the template in the wrong position (first or second relative to the target) as expected according to the modpipe style";

    $self->{ErrorCodes}->{"mismatch_cleavage_sequence_in_alignment_file"} = "Alignment files that did not have the annotated cleavage sequence at the expected position for the target sequence.";

    $self->{ErrorCodes}->{"unexpected_lysate_call"} = "Entries in lysate file that did not have 'A' or 'P' in expected column";

}

sub checkResidue{

    my ($self, $residueOneLetter, $positionInCleavageSeq, $cleavageSequence, $seqId, $peptideStartPosition) = @_;

    my @residues = split('', $cleavageSequence);

    my $codeInCleavageSeq = $residues[$positionInCleavageSeq - 1];
    if ($residueOneLetter ne $codeInCleavageSeq){
	$self->writeLog("ERROR: cleavage seq error:  cleavage seq has residue $codeInCleavageSeq at position $positionInCleavageSeq but in dssp file the residue is $residueOneLetter for sequence $seqId position $peptideStartPosition ");
    }
}



sub isError{

    my ($self, $modbaseSeqId, $cleavageSeqId, $errorCode, $modelId) = @_;
    
    my $errorDescription = $self->{ErrorCodes}->{$errorCode};
    #check to see if this combination of $modbaseSeqId and $cleavageSeqId have been reported as errors with this $errorCode
    #uses $modelId if present in the error structure too

    #check valid error code
    unless ($errorDescription){
	
	$self->writeLog("Bailing out of pipeline due to code error");

	$self->writeAllErrors();
	$self->writeLog("Code Error:  attempted to run isError() with $errorCode which is not a valid error code");
	
	exit(1);
    }
    
    my $isError = 0;
    #slightly complicated way of going through {Errors} but if we do it all at once then it populates $modbaseSeqId as key to {Errors}
    my $testModbaseSeqId =  $self->{Errors}->{$errorCode}->{$modbaseSeqId};
    if ($testModbaseSeqId){   
	my $testCleavageSeqId =  $self->{Errors}->{$errorCode}->{$modbaseSeqId}->{$cleavageSeqId};
	if ($testCleavageSeqId){
	    if ($modelId){  #if $modelId is part of the error, check to see if it is present in the structure
		my $testModel = $self->{Errors}->{$errorCode}->{$modbaseSeqId}->{$cleavageSeqId}->{modelId};
		if ($testModel){
		    $isError = 1;
		}  #else leave $isError as 0
	    }
	    else {  #no model id, automatically is an error
		$isError = 1;
	    }
	}
    }
    return $isError;
}

sub writeError{

    my ($self, $modbaseSeqId, $cleavageSeqId, $modelId, $errorCode, $methodName, $errorMessage, $writeLog) = @_;

    my $errorDescription = $self->{ErrorCodes}->{$errorCode};
    unless ($errorDescription){
	
	$self->writeLog("Bailing out of pipeline due to code error");

	$self->writeAllErrors();
	$self->writeLog("Code Error:  attempted to log error code $errorCode which is not a valid error code");
	
	exit(1);
    }
	

    $self->{Errors}->{$errorCode}->{$modbaseSeqId}->{$cleavageSeqId}->{modelId} = $modelId;
    $self->{Errors}->{$errorCode}->{$modbaseSeqId}->{$cleavageSeqId}->{errorMessage} = $errorMessage;
    $self->{Errors}->{$errorCode}->{$modbaseSeqId}->{$cleavageSeqId}->{methodName} = $methodName;
    push (@{$self->{AllPeptides}->{peptides}->{$cleavageSeqId}->{errors}}, $errorCode);

    if ($writeLog){
	$self->writeLog("ERROR:  $methodName modbase seq id: $modbaseSeqId cleavage sequence: $cleavageSeqId model id: $modelId:\n\t$errorMessage\n");
    }

}


sub writeLog {

    my ($self, $message) = @_;
    my $logFh = $self->{Log};
    my ($sec,$min,$hour,$mday,$mon,$year,$wday, $yday,$isdst)=localtime(time);

    my $dateLine = sprintf ("%4d-%02d-%02d %02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);

    print $logFh "$dateLine (Pipeline):\t$message\n";
}

sub makeTimeStamp{
    my ($self) = @_;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday, $yday,$isdst)=localtime(time);

    $year += 1900;
    $mon+=1;
    if ($mon =~ /^\d$/){
	$mon = "0" . $mon;
    }
    if ($mday =~ /^\d$/){
	$mday = "0" . $mday;
    }
    my $timeStamp = $year . "_" . $mon . "_" . $mday;
    return $timeStamp;
}


sub getRunInfo{
    my ($self) = @_;
    return $self->{RunInfo};
}



sub loadRunInfo{

    my ($self) = @_;

    my $runInfo;
    $runInfo->{"human"}->{path} = "/park1/modbase/projects/genomes/human/data/";
    $runInfo->{"human"}->{offset} = 1;
    $runInfo->{"human"}->{style} = "new";
  
    $runInfo->{"human_4-2007"}->{path} = "/park1/modbase/projects/rachel/human_4-2007/data/";
    $runInfo->{"human_4-2007"}->{offset} = 0;
    $runInfo->{"human_4-2007"}->{style} = "new";    

#    $runInfo->{"human_2008"}->{path} = "/park2/modbase/projects/genomes/human_2008/data/";
    $runInfo->{"human_2008"}->{path} = "/netapp/sali/peptide/data/landing/alignments/";

    $runInfo->{"human_2008"}->{offset} = 0;
    $runInfo->{"human_2008"}->{style} = "new";    

    $runInfo->{"snp-human2"}->{path} = "/park2/modbase/snp-human2/data/";
    $runInfo->{"snp-human2"}->{offset} = 0;
    $runInfo->{"snp-human2"}->{style} = "old";

    $runInfo->{"snp-human3"}->{path} = "/park2/rachelk/alto3/ModPipe+/snp-human3/data/";
    $runInfo->{"snp-human3"}->{offset} = 0;
    $runInfo->{"snp-human3"}->{style} = "old";

    $runInfo->{"snp-human4"}->{path} = "/park2/rachelk/netapp/rachelk/ModPipe+/snp-human4/data/";
    $runInfo->{"snp-human4"}->{offset} = 0;
    $runInfo->{"snp-human4"}->{style} = "old";

    $runInfo->{"trembl2004"}->{path} = "/park2/modbase/TrEMBL2004/data/";
    $runInfo->{"trembl2004"}->{offset} = 0;
    $runInfo->{"trembl2004"}->{style} = "old";

#    $runInfo->{"ucla-smsl"}->{path} = "";
#    $runInfo->{"ucla-smsl"}->{offset} = 0;

#    $runInfo->{"ucla-parotid"}->{path} = "";
#    $runInfo->{"ucla-parotid"}->{offset} = 0;

    $self->{RunInfo} = $runInfo;
}



sub fulfillsMethodDependencies{

    my ($self, $methodName) = @_;
     my $dependencies = $self->{MethodDependencies}->{$methodName};
    foreach my $dependency (keys %$dependencies){

	my $dependencyIsSet = $self->getParam($dependency);
	if ($dependencyIsSet eq "no"){
	    return 0;
	}
	
    }
    return 1;
}

sub writeMethodDependencies{

    my ($self) = @_;
    my $allDependencies = $self->{MethodDependencies};
    foreach my $method (keys %$allDependencies){
	my $dependencies = $allDependencies->{$method};
	my $methodString = "$method: ";
	foreach my $dependency (keys %$dependencies){
	    $methodString .= "$dependency, ";
	}
	$self->writeLog($methodString);
    }

}



sub runMethod{

    my ($self, $methodName, $skipOutput) = @_;
    
    my $runMethodParamValue = $self->getParam($methodName);

    if ($runMethodParamValue eq "yes"){
	
	if ($self->fulfillsMethodDependencies($methodName) == 1){
	    
	    $self->writeLog("Running $methodName") unless $skipOutput;
	    return 1;
	}
	else {
	    $self->writeLog("Exiting.  Pipeline step $methodName was run without having other methods run which it depended on.  Full list of method dependencies:");
	    $self->writeMethodDependencies();
	    exit(1);
	}
		
    }
    elsif ($runMethodParamValue eq "no"){
	$self->writeLog("Skipping $methodName") unless $skipOutput;
	return 0;
    }
    else {
	$self->writeLog("Please set parameter value for running method '$methodName' to either 'yes' or 'no' (case sensitive)");
	exit(1);
    }
}

sub loadMethodDependencies{
    my ($self) = @_;

    $self->{MethodDependencies}->{"runDssp"}->{"getBestModels"} = 1; 
    $self->{MethodDependencies}->{"parseDsspResults"}->{"getBestModels"} = 1; 
    $self->{MethodDependencies}->{"getTemplateLoopAlignment"}->{"getBestModels"} = 1; 

    $self->{MethodDependencies}->{"parseDsspResults"}->{"runDssp"} = 1; 

    $self->{MethodDependencies}->{"getBestModels"}->{"getSequenceLengths"} = 1;
    
#    $self->{MethodDependencies}->{"addLysateExpression"}->{"mapAccessions"} = 1;   only include if we are using secondary accessions in lysate info (currently only primary)
    $self->{MethodDependencies}->{"addAntibodyAvailability"}->{"mapSecondaryAccessions"} = 1;   

    #perhaps processNameAnnotation
    
}


sub handleReaderError{
    my ($self, $inputFile) = @_;
    my $inputFh = FileHandle->new("<" . $inputFile);
    unless ($inputFh){
	print STDERR "BenchmarkerPipeline.pl: Exiting due to FastaReader error; could not open fasta file $inputFile for reading: $!\n";
	exit(1);
    }
    else {
	print STDERR "BenchmarkerPipeline.pl: Exiting due to FastaReader error; Reader did not return true value when reading $inputFile\n";
	#this should never happen, extreme sanity check
	exit(1);
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

