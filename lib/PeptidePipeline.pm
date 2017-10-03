package PeptidePipeline;


use FastaReader;

use strict;
use FileHandle;
use POSIX qw(ceil floor);
use File::Temp;
use Cwd;

our @ISA = qw(Pipeline);

sub new{
    my ($class, $parameterFileName) = @_;
    my $self = $class->SUPER::new($parameterFileName, "peptide");
    bless $self, $class;

    $self->initialize();
    return $self;
}


sub initialize{
    my ($self) = @_;

    #error stuff
    $self->loadErrorCodes();

    $self->loadRunInfo();
}



############################################################################################################################################################                                     
# parsePeptidesFromSequences
# Top level method intended to be called when server is running in Application Scan mode. This reads the input fasta file and the rules file, and calls auxilary methods to 
# create hashes containing all sequence and peptide information. This is stored in $self->{AllSequences} which is referred to by all other methods to get
# sequences and peptides for processing.
# Format of $self->{AllSequences}:
# $allSequences->{$modbaseSeqId}->{uniprotAccession} = $uniprot_accession
#                               ->{fullSequence}     = $full_protein_sequence
#                               ->{sequenceLength}   = $full_protein_sequence_length
#                               ->{peptides}->{$peptide_start_position}->{peptideSequence}      = $peptide_sequence
#                                                                      ->{peptideStartPosition} = $peptide_start_position (this duplicates the key but makes it easier to get values later)
#                                                                      ->{peptideEndPosition}   = $peptide_end_position
#                                                                      ->{classification}       = "Application"  (other methods that load peptides have different classifications)                                             
# RETURN NULL
###########################################################################################################################################################       
sub parsePeptidesFromSequences{

    my ($self) = @_;

    $self->writeLog("Pipeline Application Mode; parsing peptides from protein sequences using rules file");
    my $startTime = time();
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $run = $self->getParam("run_name");

    #process rules file
    my $peptideRulesFileName = $self->getParam("rules_file_name");

    my $fullRulesFileName = "$pipelineDir/$run/$peptideRulesFileName";
    my $peptideLength = $self->getParam("peptide_length");
    my $peptideRules = $self->readPeptideRulesFile($fullRulesFileName);
    return 0 if $peptideRules == -1;
    $self->{PeptideRules} = $peptideRules;

    #read fasta file
    my $peptideFastaFileName = $self->getParam("input_fasta_file_name");
    my $fullFastaFileName = "$pipelineDir/$run/$peptideFastaFileName";
    my $allSequences = $self->scanPeptideFastaFile($fullFastaFileName, $peptideLength, $peptideRules);
    return 0 unless $allSequences;

    $self->{AllSequences} = $allSequences;
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "parsePeptidesFromSequence");
}


###########################################################################################################################################################       
# readPeptideInputFile
# Top level method intended to be called when server is running in Application Defined mode or Training mode. This reads the input fasta file and parses
# the headers to get the specified peptides. This is stored in $self->{AllSequences} which is referred to by all other methods to get sequences and peptides
# for processing.
#
# Format of $self->{AllSequences} is specified in parsePeptidesFromSequences() in this module.
# Format of file is Fasta, with one header line and protein sequence per accession.                                                                                                                  
# Each header line is of the form: >ModbaseSeqId|UniprotAccession(|peptideStartPosition_peptideSequence_peptideClassification)+                                                                      
# (one triplet of peptide info for each peptide, separated by '|')        
#
# RETURN NULL
###########################################################################################################################################################
sub readPeptideInputFile{

    my ($self) = @_;

    $self->writeLog("Pipeline Training Mode; reading input file to retrieve user-specified peptides");
    my $startTime = time();
    
    #load sequences using FastaReader
    my $peptideInputFileName = $self->getParam("input_fasta_file_name");
    my $pipelineDir = $self->getParam("cluster_pipeline_directory");
    my $run = $self->getParam("run_name");
    my $fullPeptideFileName = "$pipelineDir/$run/$peptideInputFileName";
    my $fastaReader = FastaReader->new($fullPeptideFileName, 1);

    unless ($fastaReader =~ /FastaReader/){
	$self->writeError("global", "global", "global", "file_missing", "readPeptideInputFile", "could not open peptide fasta input file $fullPeptideFileName: $fastaReader", 1, 1); #tested
	return 0;
    }
    
    $fastaReader->read();
    my $sequences = $fastaReader->getSequences();

    my $keywordMismatch = $self->getParam("keyword_peptide_sequence_mismatch");

    my $allSequences;

    my $sequenceCount = 0;
    my $peptideCount = 0;

    #read sequences
    foreach my $sequenceInfo (keys %$sequences){
	my $sequenceList = $sequences->{$sequenceInfo}; 
	                                                
	#parse sequence header
	my $sequence = $sequenceList->[0];
	my @cols = split('\|', $sequenceInfo);
	my $modbaseSeqId = $cols[0];
	
	my $uniprotAccession = $cols[1];
	my $peptides;
	my $peptideLength;
	$sequenceCount++;
	#get all peptides to train on from this modbaseSeqId
	my $hasOnePeptide = 0;
	for (my $i = 2; $i < scalar(@cols); $i++){
	    my $nextPeptideSpec = $cols[$i];
	    my ($peptideStartPosition, $peptideSequence, $classification) = split("\_", $nextPeptideSpec);
	    next if ($peptideSequence eq $keywordMismatch);
	    $hasOnePeptide = 1;
	    $peptideCount++;

	    $peptideLength = length($peptideSequence);
	    my $peptideEndPosition = $peptideStartPosition + $peptideLength - 1;

	    $peptides->{$peptideStartPosition}->{peptideSequence} = $peptideSequence;
	    $peptides->{$peptideStartPosition}->{peptideStartPosition} = $peptideStartPosition;   #this duplicates keys of $allSequences but makes it easier to get the values when needed later
	    $peptides->{$peptideStartPosition}->{peptideEndPosition} = $peptideEndPosition;
	    $peptides->{$peptideStartPosition}->{classification} = $classification;
	}

	if ($hasOnePeptide == 0){
	    #all sequences should have at least one peptide; this should have been validated in the front-end so is never expected to occur
	    $self->writeError("global", "global", "global", "no_defined_peptides", "readPeptideInputFile", "did not find any defined peptides for sequence $modbaseSeqId", 1, 1);
	    return 0;
	}

	#create final $allSequences data structure
	$allSequences->{$modbaseSeqId}->{uniprotAccession} = $uniprotAccession;
	$allSequences->{$modbaseSeqId}->{peptides} = $peptides;
	$allSequences->{$modbaseSeqId}->{fullSequence} = $sequence;
	my $sequenceLength = length($sequence);

	$allSequences->{$modbaseSeqId}->{sequenceLength} = $sequenceLength;
    }
    
    $self->writeLog("Done loading peptides; loaded $peptideCount peptides from $sequenceCount sequences");

    $self->{AllSequences} = $allSequences;
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "readPeptideInputFile");
}


###########################################################################################################################################################       
# scanPeptideFastaFile
# auxilary method that does the work for parsePeptidesFromSequence.
# 
# PARAM  $peptideFastaFile: input file in fasta format containing sequence info. See peptide.pm in saliweb frontend (method writeApplicationInfo) for exact format.
# PARAM  $peptideLength: length of the peptides to be parsed. 
# PARAM  $peptideRules: hash containing allowed residues at each position (format specified in readPeptideRulesFile())
# RETURN $allSequences: hash containing all info extracted from fasta file (format specified in parsePeptidesFromSequence())
###########################################################################################################################################################       
sub scanPeptideFastaFile{

    my ($self, $peptideFastaFile, $peptideLength, $peptideRules) = @_;

    #load sequences using FastaReader
    my $fastaReader = FastaReader->new($peptideFastaFile, 1);

    unless ($fastaReader =~ /FastaReader/){
	$self->writeError("global", "global", "global", "file_missing", "scanPeptideFastaFile", "could not open peptide fasta input file $peptideFastaFile: $fastaReader", 1, 1); 
	return 0;
    }
    $fastaReader->read();
    my $sequences = $fastaReader->getSequences();
    
    #go through headers and sequences, load into $allSequences
    my $allSequences;

    my $sequenceCount = 0;
    my $peptideCount = 0;
    foreach my $sequenceInfo (keys %$sequences){   #$sequenceInfo: header in fasta file of format $modbaseSeqId|$uniprotAccession
	my $sequenceList = $sequences->{$sequenceInfo};   
	my $sequence = $sequenceList->[0];
	my ($modbaseSeqId, $uniprotAccession) = split('\|', $sequenceInfo);
	#parse peptides using rules file
	my $peptides = $self->getPeptidesFromFullSequence($sequence, $peptideLength, $peptideRules);
	
	#load info
	$allSequences->{$modbaseSeqId}->{uniprotAccession} = $uniprotAccession;
	$allSequences->{$modbaseSeqId}->{peptides} = $peptides;
	$allSequences->{$modbaseSeqId}->{fullSequence} = $sequence;
	my $sequenceLength = length($sequence);
	$allSequences->{$modbaseSeqId}->{sequenceLength} = $sequenceLength;
	
	#stats
	my $peptideCountForSequence =  scalar(keys %$peptides);
	$peptideCount += $peptideCountForSequence;
	if ($peptideCountForSequence == 0){
	    $self->writeNoPeptidesParsedForSeq($modbaseSeqId);
	    $self->writeLog("Did not parse peptides from modbase sequence $modbaseSeqId\n");
	}
	$sequenceCount++;
    }

    $self->writeLog("Done loading peptides; loaded $peptideCount peptides from $sequenceCount sequences");

    return $allSequences;
}


###########################################################################################################################################################       
# getPeptidesFromFullSequence
# Use peptide rules to parse peptides from full protein sequence.
# PARAM  $sequence: full protein sequence input
# PARAM  $peptideLength: length of the peptides to be parsed. 
# PARAM  $peptideRules: hash containing allowed residues at each position (format specified in readPeptideRulesFile())
# RETURN $peptides: hash containing all peptide info (format specified in parsePeptidesFromSequences())
###########################################################################################################################################################       
sub getPeptidesFromFullSequence{
    
    my ($self, $sequence, $peptideLength, $peptideRules) = @_;
    
    my $peptides;
    my @seqArray = split('', $sequence);

    for (my $i = 0; $i < scalar(@seqArray) - $peptideLength + 1; $i++){  
	my $nextPeptide = substr($sequence, $i, $peptideLength);

	#add all peptide info to $peptides, unless this particular peptide should be filtered according to rules
	if ($self->passesCheck($nextPeptide, $peptideRules)){     
	    $peptides->{$i}->{peptideSequence} = $nextPeptide;
	    my $peptideEndPosition = $i + $peptideLength - 1;
	    $peptides->{$i}->{peptideStartPosition} = $i;   #this duplicates keys of $allPeptides but makes it easier to get the values when needed later
	    $peptides->{$i}->{peptideEndPosition} = $peptideEndPosition;
	    $peptides->{$i}->{classification} = "Application";
	}
    }
    return $peptides;
}

###########################################################################################################################################################       
# passesCheck
# Helper method for getPeptidesFromFullSequence; examines one peptide and determines if any positions contain a residue that the user doesn't want present
#
# PARAM  $peptide: The peptide in question
# PARAM  $peptideRules: hash containing allowed residues at each position (format specified in readPeptideRulesFile())
# RETURN 0 if the peptide has a residue that the user doesn't want present; 0 otherwise
###########################################################################################################################################################       
sub passesCheck{
    
    my ($self, $peptide, $peptideRules) = @_;
    foreach my $positionNumber (keys %$peptideRules){
	
	#read residues that cannot be present at this position
	my $bannedResidues = $peptideRules->{$positionNumber};
	my $residueAtPosition = substr($peptide, $positionNumber - 1, 1);  #position is 1-based

	#check if any of these residues are indeed present
	foreach my $bannedResidue (keys %$bannedResidues){
	    $bannedResidue = uc($bannedResidue);
	    if ($bannedResidue eq $residueAtPosition){
		return 0;
	    }
	}
    }
    return 1;
}


############################################################################################################################################################                                         
# readPeptideRulesFile
# Reads 'rules file' which specifies the rules for how to find peptides given a full protein sequence. 
# Full details of rules file format are in the peptide.pm saliweb frontend module.
#
# PARAM  $peptideRulesFileName: full path to the rules file
# RETURN hash of the form: $peptideRules->{$peptidePositionNumber}->{$oneLetterResidue} = 1
#                          $oneLetterResidue is a residue that should NOT be present in any peptide. So if you only want Asp to be present 
#                          in the fourth position of a peptide eg with GrB, all residues except for D would be keys of $peptideRules->{4}->{}
############################################################################################################################################################      
sub readPeptideRulesFile{
    my ($self,  $peptideRulesFileName) = @_;


    my $peptideRulesFh = FileHandle->new("<" . $peptideRulesFileName);
    unless ($peptideRulesFh =~ /FileHandle/){
	$self->writeError("global", "global", "global", "file_missing", "readPeptideRulesFile", "could not open peptide rules file $peptideRulesFileName for reading: $!", 1, 1); #tested
	return -1;
    }

    my $peptideRules;

    while (<$peptideRulesFh>){

	chomp;
	my $line = $_;               #line format: 1 A C D F  where 1 is the position in the peptide to restrict and is followed by list of residues
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



############################################################################################################################################################      
# getBestModels
# Get the best scoring models for each peptide. The peptide must be fully contained within the model, so we iterate through all models for each sequence
# to see which contain the peptide.  Then of those models that do contain it, choose the one that scores the best according to user-specified best model criteria.
# This criteria can either be the ModPipe model score, coverage of the model relative to the entire protein, or TSVMod score.
#
# Stores the results (i.e. all information of the best model for a peptide) in 
# $self->{AllSequences}->{$modbaseSeqId}->{$peptideStartPosition}->{$bestModelId}->{info_name} = $infoValue
# <infoName> keys are a dozen or so model info values (start and end residues of model coverage, dope score, etc.). 
#
# Also stores other information about the modeling process (templates, model count, etc.) in 
# $self->{AllSequences}->{$modbaseSeqId}->{$peptideStartPosition}
#############################################################################################################################################################      
sub getBestModels{

    my ($self) = @_;
    
    $self->writeLog("Getting best models for sequences");

    return 0 if ($self->pipelineHasErrors());
	
    my $startTime = time();

    my $modelTable = $self->getParam('model_table');

    my $totalModelsForSequencesCount = 0;
    my $totalModelsContainingPeptidesCount = 0;
    my $bestModelCount = 0;
    my $sequenceCount = 0;
    my $uniqueSeqsWithModels;

    my $allSequences = $self->{AllSequences};

    my $bestModelCriteria = $self->getParam("best_model_criteria");    
    foreach my $modbaseSeqId (keys %$allSequences){

	my $fh = FileHandle->new("<" . $modelTable);
	
	unless ($fh =~ /FileHandle/){
	    $self->writeError("global", "global", "N/A", "file_missing", "getBestModels", "could not open model table $modelTable: $!", 1, 1); #tested
	    return 0;
	}

	$sequenceCount++;

	my $peptides = $allSequences->{$modbaseSeqId}->{peptides};
	my $sequenceLength = $allSequences->{$modbaseSeqId}->{sequenceLength};
	my $peptideLength = $self->getParam("peptide_length");
	
	my $modelInfo;
	
	#read through model file, get info for all models for this sequence
	my $modelsForSequenceCount = 0;
	my $modelErrors;
	while (<$fh>){
	    chomp;
	    my $line = $_;
	    my @cols = split('\t', $line);
	    my $seqId = $cols[0];

	    if ($seqId eq $modbaseSeqId){

		$modelsForSequenceCount++;
		$totalModelsForSequencesCount++;
		#read model input file (results of "select * from newmodels where run = 'human_2008'")
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

		#calculate coverage for model
		my $modelLength = $modelTargetEnd - $modelTargetStart;
		if ($modelLength < 1){
		    my $msg = "Invalid model table $modelTable (modbase seq $modbaseSeqId model id $modelId). Got invalid model length of $modelLength (model start: $modelTargetStart, model end: $modelTargetEnd)";
		    $self->writeError("global", "global", $modelId, "invalid_model", "getBestModels", $msg, 1, 1); #tested
		    last;
		}
		my $coverage = ($modelLength * 1.0) / ($sequenceLength * 1.0);
		my $roundedCoverage = sprintf("%.2f", $coverage);
		
		#store all model info for this id
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
	$fh->close();
	$allSequences->{$modbaseSeqId}->{totalModelCount} = $modelsForSequenceCount;
	
	#for each peptide, find best model that contains it
	foreach my $peptideStartPosition (keys %$peptides){
	    my $peptideInfo = $peptides->{$peptideStartPosition};
	    my $peptideEndPosition = $peptideInfo->{peptideEndPosition};
	    my $modelsContainingThisPeptideCount;
	    my @modelsContainingPeptide;

	    #add all models containing peptide to @modelsContainingPeptide
	    foreach my $modelId (keys %$modelInfo){
		my $modelTargetStart = $modelInfo->{$modelId}->{modelTargetStart};
		my $modelTargetEnd = $modelInfo->{$modelId}->{modelTargetEnd};
		
		if ($modelTargetStart <= $peptideStartPosition && $modelTargetEnd >= $peptideEndPosition){
		    push(@modelsContainingPeptide, $modelId);
		    $modelsContainingThisPeptideCount++;
		    $totalModelsContainingPeptidesCount++;
		}
	    }
	    
	    #use model properties and critiera for best model (user specified) to get the model id of the highest quality model for this peptide
	    my $bestModelId = $self->getBestModelByCriteria($bestModelCriteria, $modelInfo, @modelsContainingPeptide);
	    my $modelTargetStart = $modelInfo->{$bestModelId}->{modelTargetStart};
	    my $modelTargetEnd = $modelInfo->{$bestModelId}->{modelTargetEnd};
	    
	    #put all info for best model into $peptideInfo
	    if ($bestModelId){
		
		foreach my $infoType (keys %{$modelInfo->{$bestModelId}}){
		    $peptideInfo->{$infoType} = $modelInfo->{$bestModelId}->{$infoType};
		}
		$peptideInfo->{bestModelId} = $bestModelId;
		$bestModelCount++;
		$uniqueSeqsWithModels->{$modbaseSeqId} = 1;
		my $alignmentId = $peptideInfo->{alignmentId};
		my $run = $self->getParam("modpipe_run_name");
		
		my $fullAlignmentFileName = $self->getAlignmentFileName($modbaseSeqId, $alignmentId, $run);
		
		my $templatePdbIds = $self->getTemplatePdbId($modbaseSeqId, $peptideStartPosition, $bestModelId, $fullAlignmentFileName, $run);
		
		$peptideInfo->{templatePdbIds} = $templatePdbIds;
		$peptideInfo->{fullAlignmentFilePath} = $fullAlignmentFileName;
		$peptideInfo->{modelsContainingPeptideCount} = $modelsContainingThisPeptideCount; 
	    }
	}
    }

    my $uniqueSeqCount = scalar (keys %$uniqueSeqsWithModels);

    my $outputString ="Done getting best models. Processed $sequenceCount sequences, found $totalModelsForSequencesCount total models in $uniqueSeqCount unique sequences, $totalModelsContainingPeptidesCount models containing at least one peptide, "; 
    $outputString .= " and $bestModelCount peptides had a best-scoring model according to criteria $bestModelCriteria";
    $self->writeLog($outputString);

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "getBestModels");
}


############################################################################################################################################################      
# getAlignmentFileName
# Get the full path and file name of the alignment file for the given alignment ID. Does not check as to whether file exists (this is done by other methods that
# attempt to read the file).
# 
# PARAM  $modbaseSeqId: modbase sequence id for the alignment
# PARAM  $alignmentId:  id for the alignment (checksum format)
# PARAM  $runName:      name of the modpipe run that was used to generate the alignment
# RETURN $alignmentFileName: String with full path of file name
############################################################################################################################################################      
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


############################################################################################################################################################      
# getTemplatePdbId
# Read given alignment file and get the PDB codes for the templates that were aligned against the given sequence.
#
# PARAM  $modbaseSeqId: sequence id for the target in the alignment
# PARAM  $peptideStartPosition: starting position for the peptide being processed (only for error output purposes)
# PARAM  $fullAlignmentFileName: full alignment file path and name from which templates will be read
# PARAM  $run: name of the modpipe run that generated the alignment
# RETURN $pdbCode: 4-letter PDB code for the template
############################################################################################################################################################      
sub getTemplatePdbId{

    my ($self, $modbaseSeqId, $peptideStartPosition, $modelId, $fullAlignmentFileName, $run) = @_;
    
    my $runInfo = $self->getRunInfo();
    my $runStyle = $runInfo->{$run}->{style};

    #read alignment file to get template pdb ids
    my ($targetInfoLine, $targetSequence, $templateInfoLine, $templateSequence) = $self->parseAlignmentFile($fullAlignmentFileName, $run, $modbaseSeqId, $peptideStartPosition, $modelId);
    
    return "" unless $targetInfoLine;

    my @cols = split('\:', $templateInfoLine);
    my $pdbCode = $cols[1];
    my $proteinType = $cols[0];  #either 'structure' or 'structureX' depending on run type

    #check consistency in modpipe alignment file format
    unless ($proteinType =~ /structure/){   
	my $errorString = "Did not get protein type of 'structure' (got $proteinType) in first column of pir format in expected position in alignment file $fullAlignmentFileName run style $runStyle";
	$self->writeError($modbaseSeqId, $peptideStartPosition, 0, "incorrect_template_position_in_alignment", "getBestModels", $errorString, 1, 0); #Untested (after $runInfo parameterized)
	return "";
    }
   
    unless ($pdbCode =~ /^\S\S\S\S$/){    #expect 4 letter pdb code
	my $errorString = "Did not get four letter pdb code (got $pdbCode) in first column of pir line in expected position in alignment file $fullAlignmentFileName run style $runStyle";
	$self->writeError($modbaseSeqId, $peptideStartPosition, 0, "unexpected_pdb_code_in_alignment", "getBestModels", $errorString, 1, 0); #Untested (after $runInfo parameterized)
	return "";
    }

    return $pdbCode;
}



############################################################################################################################################################      
# getBestModelByCriteria
# Reads the list of models and returns the one that scores the highest according to the specified criteria.
# A bit of a hack since all criteria assume the best is always the one with the greatest score; this is not true for scores such as DOPE (which currently aren't
# supported best criteria, but if included, would necessitate changing the flow). In case of a tie, returns the first model in the list with the tied score.
# 
# PARAM  $bestModelCriteria: string from a CV of possible criteria for best model. One of "nativeOverlap", "coverage", "modelScore"
# PARAM  $modelInfo: hash where the keys are $modelIds and values are all information about the model, including the one to find the given criteria
# PARAM  @modelsContainingPeptide: array of model IDs under consideration
# RETURN model id (checksum form) that has the best score
############################################################################################################################################################      
sub getBestModelByCriteria{

    my ($self, $bestModelCriteria, $modelInfo, @modelsContainingPeptide) = @_;   

    my $currentBestCriteria = 0.0; #assumes that the minimum score for all criteria is >= 0
    my $currentBestModelId;
    foreach my $modelId (@modelsContainingPeptide){
	my $currentCriteria = $modelInfo->{$modelId}->{$bestModelCriteria};  #this is assumed to have been populated already (done in getBestModels() which calls this method)
	if ($currentCriteria > $currentBestCriteria){
	    $currentBestModelId = $modelId;
	    $currentBestCriteria = $currentCriteria;
	}
    }
    return $currentBestModelId;
}


############################################################################################################################################################      
# parseDsspResults 
# For all peptides that have models associated with them, read DSSP file to get structure elements associated with the peptides in the models
# DSSP files have been pre-calculated for all human models in 2008 run (will need to be updated accordingly). Results go into
# $self->{AllSequences}->{$modbaseSeqId}->{$peptideStartPosition}->{featureType} = $featureValue
# (featureType details in method comments)
#
# RETURN Null
############################################################################################################################################################ 

sub parseDsspResults{
    my ($self) = @_;
    
    $self->writeLog("Parsing results of DSSP runs to get structure information from structures and models");
    return 0 if ($self->pipelineHasErrors());
    my $startTime = time();

    my $dsspDir = $self->getParam("dssp_directory");  

    my $allSequences = $self->{AllSequences};
    my $residueSolventAcc = $self->getResidueSolventAcc();
    
    my $sequenceCount = 0;
    my $peptideCount = 0;
    my $peptidesWithDsspCount = 0;
    my $errorCount = 0;

    foreach my $modbaseSeqId (keys %$allSequences){
	my $peptides = $allSequences->{$modbaseSeqId}->{peptides};
	$sequenceCount++;
	
	foreach my $peptideStartPosition (keys %$peptides){
	    
	    my $peptideInfo = $peptides->{$peptideStartPosition};
	    $peptideCount++;
	    my $peptideSequence = $peptideInfo->{peptideSequence};
	    
	    #see if this peptide has a model that contains it (set in getBestModels())
	    $modbaseSeqId =~ /^(\S\S\S)\S+/;
	    my $seqTopDir = $1;
	    my $bestModelId = $peptideInfo->{bestModelId};
	    next unless $bestModelId;
	
	    #get dssp file for this model (pre-calculated)
	    my $noFile = 0;
	    my $dsspFile = "$dsspDir/$seqTopDir/$modbaseSeqId/$bestModelId.dssp";
	    my $dsspFh = FileHandle->new("<" . $dsspFile);
	    unless ($dsspFh =~ /FileHandle/){
		$self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId, "no_dssp_structure_info", "parseDsspResults", "Could not open dssp file $dsspFile: $!", 1, 0); #tested
		$errorCount++;
		next;
	    }
	    
	    #Model residue numbering varies according to modpipe run style. Account for that here to get peptide range as it appears in DSSP file
	    my ($modelPeptideStartPosition, $modelPeptideEndPosition) = $self->getModelPeptideRange($peptideInfo);
	    if ($modelPeptideStartPosition =~ /\D/){ #error message
		my $msg = $modelPeptideStartPosition;
		$self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId, "modpipe_run_info_error", "parseDsspResults",  $msg, 1, 0);
		$errorCount++; ##Untested (after $runInfo parameterized)
		next;
	    }
	    #convert from 0 to 1 based (DSSP file is 1-based)   
	    $modelPeptideStartPosition++;
	    $modelPeptideEndPosition++;     

	    #read structure info
	    my $peptideSiteSSInfo;
	    my $currentResiduePosition = 1; #starts at one no matter which style modpipe was used
	    my $baseOnePeptideStart;  #peptide sites if pdb file were counted starting at one
	    my $baseOnePeptideEnd;
	    my @structureValueArray;
	    my @structureTypeArray;
	    my @accessibilityFractionArray;
	    my @accessibilityTypeArray;
	    my $skipProcessing = 0;	    
	    
	    #move down past header info
	    while (<$dsspFh>){
		chomp;
		my $line = $_;
		last if ($line =~ /\s+\#/);
	    }

	    while (<$dsspFh>){
		chomp;
		my $line = $_;
		
		last if ($line =~ /^\s+$/);
		
		# Format of line: 1135 1145   G  S >  S+     0   0   56  .*
		# Description of each column at http://swift.cmbi.kun.nl/gv/dssp/		                   
		if ($line =~ /\s+\d+\s+(\d+)\s+(\w)\s\s(.).{12}\s*\S+\s+(\d+).*/){
		    
		    my $dsspResidueNumber = $1;
		    my $residueOneLetter = $2;
		    my $structureType = $3;
		    my $accessibility = $4;
		    $currentResiduePosition++;
		    
		    $peptideSiteSSInfo->{$currentResiduePosition} = $structureType;  #prepare for counting loop length
		    
		    #check if current residue is within the peptide sequence boundaries
		    if ($dsspResidueNumber >= $modelPeptideStartPosition && $dsspResidueNumber <= $modelPeptideEndPosition){
			
			#make sure peptide sequence in dssp file matches $peptideSequene for peptide we are processing
			my $positionInPeptideSeq = $dsspResidueNumber - $modelPeptideStartPosition + 1;
			my $match = $self->checkResidueMatchesDssp($residueOneLetter, $positionInPeptideSeq, $peptideSequence, $modbaseSeqId, $modelPeptideStartPosition);

			unless ($match == 1){
			    my $errorMsg = $match;
			    $self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId, "dssp_format_error", "parseDsspResults",  $errorMsg, 1, 0);  #tested
			    $skipProcessing = 1;
			    $errorCount++;
			    last;
			}
			
			#process secondary structure
			my ($mappedStructureType, $mappedStructureValue) = $self->getMappedStructureType($structureType);
			unless ($mappedStructureType){
			    my $errorMsg = "Did not get mapped structure type for dssp structure $structureType (line $line)"; 
			    $self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId,  "dssp_format_error", "parseDsspResults", $errorMsg, 1, 0);  #tested
			    $skipProcessing = 1;
			    $errorCount++;
			    last;
			}
			push (@structureTypeArray, $mappedStructureType);
			push (@structureValueArray, $mappedStructureValue);
			
			#process solvent accessibility
			my ($percentAccessible, $accCall) = $self->getSolventExposureFraction($accessibility, $residueOneLetter, $residueSolventAcc);
			
			unless ($accCall){
			    my $errorMsg = "Did not get solvent exposure fraction for residue $residueOneLetter (line $line)";
			    $self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId,  "dssp_format_error", "parseDsspResults", $errorMsg, 1, 0); #tested
			    $skipProcessing = 1;
			    $errorCount++;
			    last;
			}
			push (@accessibilityFractionArray, $percentAccessible);
			push (@accessibilityTypeArray, $accCall);
		    }
		    
		    #prepare for counting loop length later
		    if ($dsspResidueNumber == $modelPeptideStartPosition){
			$baseOnePeptideStart = $currentResiduePosition;
		    }
		    if ($dsspResidueNumber == $modelPeptideEndPosition){
			$baseOnePeptideEnd = $currentResiduePosition;
		    }
		    
		}
		else {
		    $self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId,  "regex_error", "parseDsspResults",   #tested
				      "Did not correctly parse regex from line $line in dssp file $dsspFile", 1, 0);
		    $skipProcessing = 1;
		    $errorCount++;
		    last;
		}
	    }
	    next if $skipProcessing;

	    #create strings from structure / accessibility values
	    my $structureTypeString = join('', @structureTypeArray);
	    if ($structureTypeString eq ""){  #if we couldn't get a result for secondary structure, bail on solvent accessibility
		my $errorString = "Could not read DSSP secondary structure information from DSSP result file $dsspFile";
		$errorCount++;
		$self->writeError($modbaseSeqId, $peptideStartPosition, $bestModelId, "no_dssp_structure_info", "parseDsspResults", $errorString, 1); #tested
	    }
	    else {
	    
		my $structureValueString = join(',', @structureValueArray);
		my $accessibilityTypeString = join('', @accessibilityTypeArray);
		my $accessibilityFractionString = join(',', @accessibilityFractionArray);
		
		#set structure / accessibilty values for this peptide
		$peptideInfo->{structureTypes} = $structureTypeString;                   #structure type = "A", "L", "B" 
		$peptideInfo->{structureValues} = $structureValueString;                 #structure value = 1, 2, 3
		$peptideInfo->{accessibilityTypes} = $accessibilityTypeString;           #accessibility type = "A", "S"
		$peptideInfo->{accessibilityFractions} = $accessibilityFractionString;   #accessibility fraction = (comma separated string with fraction list)
	    
		# Previously doing this to check how long the length of the loop was if this peptide was on a loop
		#	my $loopCount = $self->countCleavageSiteLoopLength($peptideSiteSSInfo, $baseOnePeptideStart, $baseOnePeptideEnd); 
		
		my $loopCount = "Loop not counted"; 
		$peptideInfo->{"loopLength"} = $loopCount;
		$peptidesWithDsspCount++;
	    }
	}
    }

    my $logMsg = "Done loading DSSP results. Processed $sequenceCount sequences, $peptideCount peptides, and $peptidesWithDsspCount peptides found in a model for which a DSSP file was generated. ";
    $logMsg .= "$errorCount peptides were in DSSP files which were not correctly processed.";
    $self->writeLog($logMsg);

    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "parseDsspResults");
}

############################################################################################################################################################      
# getMappedStructureType
# Given the DSSP result of what element of secondary structure a residue is on, convert to either L (Loop), A (Alpha-Helix), or B (Beta-Sheet).
# DSSP has different types of helix and sheets, here we compact them into a global structure type.
#
# PARAM  $structureInputType: one letter (or space for some loops) representing DSSP's call as to the structure of the residue
# RETURN Two-element array; first entry is the letter representing the mapped call (A, L, or B) and the second is the numeric value for the call that will
#        serve as input to the SVM that processes this.
############################################################################################################################################################      
sub getMappedStructureType{

    my ($self, $structureInputType) = @_;

    if ($structureInputType =~ /^\s+$/){  #space is considered a loop or unassigned. Could also be parsing error (would be caught by checkResidueMatchesDssp)
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
    #error if we get here (handle in parseDsspResults())
    return 0;
}



############################################################################################################################################################      
# runPsipred
# For all peptides, get the Psipred call for their predicted secondary structure. If there is no Psipred result file, run it on the fly, and save results
# in data directory for later runs. Results go into 
# $self->{AllSequences}->{$modbaseSeqId}->{$peptideStartPosition}->{psipredString} = $psipredStructureType (single string)
#                                                                ->{psipredFirstScores} = $psipredScores   (string with comma separated list of scores)
#
# RETURN Null
############################################################################################################################################################      
sub runPsipred{

    my ($self) = @_; 

    my $startTime = time();
    $self->writeLog("Running Psipred to get sequence-based secondary structure information");
    return 0 if ($self->pipelineHasErrors());

    my $allSequences = $self->{AllSequences};
    my $psipredResultsDir = $self->getParam("psipred_results_dir");
    my $sequenceCount = 0;
    my $peptideCount = 0;
    my $errorCount = 0;

    my $psipredStructureMap = $self->makePsipredStructureMap();

    foreach my $modbaseSeqId (keys %$allSequences){
	$sequenceCount++;
	my $peptides = $allSequences->{$modbaseSeqId}->{peptides};
	foreach my $peptideStartPosition (keys %$peptides){

	    my $peptideInfo = $peptides->{$peptideStartPosition};

	    my $peptideStartPosition = $peptideInfo->{peptideStartPosition};
	    my $psipredStartPosition = $peptideStartPosition + 1;     #psipred expects 1-based residue numbering and ours go in as 0-based.   
	    
	    my $peptideSequence = $peptideInfo->{peptideSequence};
	    my @peptideArray = split('', $peptideSequence);
    	    my $peptideLength = scalar(@peptideArray);
	    #Get psipred results file from modbase sequence id and input parameters
	    $modbaseSeqId =~ /^(\S\S)\S+/;
	    my $twoLetters = $1;
	    my $fullPsipredDir = $psipredResultsDir . "/" . $twoLetters;
	    my $fullPsipredFile = $fullPsipredDir . "/" . $modbaseSeqId . ".ss2";

	    my $psipredResult;

	    #check if file exists; if not, run psipred
	    if (-e $fullPsipredFile == 0){
		$self->writeLog("Did not find Psipred file for modbase sequence ID $modbaseSeqId; running psipred command to create it");
		$psipredResult = $self->runPsipredCommand($modbaseSeqId, $peptideStartPosition, $fullPsipredDir);
		next if ($psipredResult == -1);
	    }
	    
	    #open psipred result file
	    my $psipredFh = FileHandle->new("<" . $fullPsipredFile);
	    unless ($psipredFh =~ /FileHandle/){
		my $msg = "Could not open psipred output file $fullPsipredFile: $!";
		$msg .= " This followed a psipred run to generate new psipred files; output from there: $psipredResult" if $psipredResult;
		$self->writeError($modbaseSeqId, $peptideStartPosition, "N/A",  "no_psipred_result_file", "runPsipred", $msg, 1, 0);  #tested
		$errorCount++;
		next;
	    }
	    
	    
	    my $structureString = "";
	    my @firstScoreArray;
	    
	    #scroll down to start of results
	    for (my $i = 0; $i < 2; $i++){
		<$psipredFh>;
	    }
	    
	    my $skipProcessing = 0;
	    while (<$psipredFh>){
		chomp;
		my $line = $_;
		if ($line =~ /^\s+(\S.*)/){
		    $line = $1;   #trim initial blank space from line for three digit residue numbers (four digit not an issue)
		}
		my ($residueNumber, $residueOneLetter, $structureCall, $firstScore, $secondScore, $thirdScore)  = split('\s+', $line);  
		
		if ($residueNumber == $psipredStartPosition){  
		    #found peptide start position, process the whole peptide

		    for (my $i = 0; $i < $peptideLength; $i++){

			#make sure peptide sequence in psipred file matches $peptideSequence for the peptide we are processing
			my $expectedPeptideResidue = $peptideArray[$i];
			if ($expectedPeptideResidue ne $residueOneLetter){
			    my $errorString =  "Error: residue $residueOneLetter in psipred results file (position $residueNumber) does not equal expected residue entry $expectedPeptideResidue";
			    $self->writeError($modbaseSeqId, $peptideStartPosition, "N/A", "psipred_sequence_mismatch", "processPsipredResults", $errorString, 1, 0); #tested
			    $skipProcessing = 1;
			    $errorCount++;
			    last;
			}
		
			#get structure type for this psipred call
			my $nextStructureType = $psipredStructureMap->{$structureCall};
			unless ($nextStructureType){
			    my $errorString = "Error: did not get mapped structure type for psipred call $structureCall (line $line)";
			    $self->writeError($modbaseSeqId, $peptideStartPosition, "N/A", "psipred_format_error", "processPsipredResults", $errorString, 1, 0); #tested
			    $skipProcessing = 1;
			    $errorCount++;
			    last;
			}

			$structureString .= $nextStructureType;
			push(@firstScoreArray, $firstScore);

			#read in next line
			$line = <$psipredFh>;
			chomp $line;
			if ($line =~ /^\s+(\S.*)/){
			    $line = $1;   #trim initial blank space from line for three digit residue numbers (four digit not an issue)
			}
			($residueNumber, $residueOneLetter, $structureCall, $firstScore, $secondScore, $thirdScore)  = split('\s+', $line); 
		    } #End peptide processing
		    
		    last; #we have finished the peptide.  Theoretically if the protein has more than one peptide we could keep going but will just process this sequence again instead
		}
	    }
	    
	    next if $skipProcessing; #proceed to next peptide
	    
	    #add results to $peptideInfo
	    if ($structureString eq ""){   #did not find expected peptide start position in psipred result file
		$self->writeError($modbaseSeqId, $peptideStartPosition, "N/A", "no_psipred_peptide_sequence", "processPsipredResults", 
				  "Found Psipred results file but did not find peptide sequence start position $peptideStartPosition in file", 1, 0); #tested
		$errorCount++;
	    }
	    else {
		$peptideCount++;
		$peptideInfo->{psipredString} = $structureString;
		my $firstScoreString = join (',', @firstScoreArray);
		$peptideInfo->{psipredFirstScores} = $firstScoreString;
	    }
	}
    }

    $self->writeLog("Done running PsiPred. Processed $sequenceCount sequences containing $peptideCount peptides. Encountered $errorCount errors during this run");
    
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "runPsipred");
}


############################################################################################################################################################      
# runDisopred
# For all peptides, get the Disopred call for their predicted disorder type. If there is no Disopred result file, run it on the fly, and save results
# in data directory for later runs. Results go into 
# $self->{AllSequences}->{$modbaseSeqId}->{$peptideStartPosition}->{disopredString} = $disopredDisorderType (single string)
#                                                                ->{disopredFirstScores} = $disopredScores   (string with comma separated list of scores)
#
# RETURN Null
############################################################################################################################################################      
sub runDisopred{

    my ($self) = @_; 

    my $startTime = time();
    $self->writeLog("Running Disopred to get sequence-based disorder predictions");
    return 0 if ($self->pipelineHasErrors());

    my $allSequences = $self->{AllSequences};
    my $disopredResultsDir = $self->getParam("disopred_results_dir");
    my $sequenceCount = 0;
    my $peptideCount = 0;
    my $errorCount = 0;

    foreach my $modbaseSeqId (keys %$allSequences){
	$sequenceCount++;
	my $peptides = $allSequences->{$modbaseSeqId}->{peptides};
	foreach my $peptideStartPosition (keys %$peptides){

	    my $peptideInfo = $peptides->{$peptideStartPosition};
	    my $peptideStartPosition = $peptideInfo->{peptideStartPosition};
	    my $disopredStartPosition = $peptideStartPosition + 1;     #disopred expects 1-based residue numbering and ours go in as 0-based.   
	    
	    my $peptideSequence = $peptideInfo->{peptideSequence};
	    my @peptideArray = split('', $peptideSequence);
    	    my $peptideLength = scalar(@peptideArray);	    
	    	    	    
	    #Get disopred results file from modbase sequence id and input parameters
	    $modbaseSeqId =~ /^(\S\S)\S+/;
	    my $twoLetters = $1;
	    my $fullDisopredDir = $disopredResultsDir . "/" . $twoLetters;
	    my $fullDisopredFile = $fullDisopredDir . "/" . $modbaseSeqId . ".diso";

	    my $disopredResult;

	    #check if file exists; if not, run disopred
	    if (-e $fullDisopredFile == 0){
		$self->writeLog("Did not find Disopred file for modbase sequence ID $modbaseSeqId; running disopred command to create it"); 
		$disopredResult = $self->runDisopredCommand($modbaseSeqId, $peptideStartPosition, $fullDisopredDir);
		next if ($disopredResult == -1);
	    }
	    
	    #open disopred result file
	    my $disopredFh = FileHandle->new("<" . $fullDisopredFile);
	    unless ($disopredFh =~ /FileHandle/){
		my $msg = "Could not open disopred output file $fullDisopredFile: $!";
		$msg .= " This followed a disopred run to generate new disopred files; output from there: $disopredResult" if $disopredResult;
		$self->writeError($modbaseSeqId, $peptideStartPosition, "N/A",  "no_disopred_result_file","runDisopred", $msg, 1, 0); #tested
		$errorCount++;
		next;
	    }
	    
	    
	    my $disorderString = "";
	    my @firstScoreArray;
	    
	    #scroll down to start of results
	    for (my $i = 0; $i < 5; $i++){
		<$disopredFh>;
	    }
	    
	    my $skipProcessing = 0;
	    while (<$disopredFh>){
		chomp;
		my $line = $_;
		my ($blank, $residueNumber, $residueOneLetter, $disorderCall, $firstScore, $secondScore)  = split('\s+', $line);  
		
		if ($residueNumber == $disopredStartPosition){  
		    #found peptide start position, process the whole peptide

		    for (my $i = 0; $i < $peptideLength; $i++){

			#make sure peptide sequence in disopred file matches $peptideSequence for the peptide we are processing
			my $expectedPeptideResidue = $peptideArray[$i];
			if ($expectedPeptideResidue ne $residueOneLetter){
			    my $errorString =  "Error: residue $residueOneLetter in disopred results file (position $residueNumber) does not equal expected residue entry $expectedPeptideResidue";
			    
			    $self->writeError($modbaseSeqId, $peptideStartPosition, "N/A", "disopred_sequence_mismatch", "processDisopredResults", $errorString, 1, 0); #tested
			    $skipProcessing = 1;
			    $errorCount++;
			    last;
			}
		
			#get disorder value for this disopred call
			if ($disorderCall eq "*"){
			    $disorderString .= "D";
			}
			elsif ($disorderCall eq "."){
			    $disorderString .= "O";
			}
			else {
			    my $errorString = "Error: did not get mapped disorder value for disopred call $disorderCall (expecting  '*' or '.'; line $line)";
			    $self->writeError($modbaseSeqId, $peptideStartPosition, "N/A", "disopred_format_error", "processDisopredResults", $errorString, 1, 0); #tested
			    $skipProcessing = 1;
			    $errorCount++;
			    last;
			}

			push(@firstScoreArray, $firstScore);

			#read in next line
			$line = <$disopredFh>;
			chomp $line;
			($blank, $residueNumber, $residueOneLetter, $disorderCall, $firstScore, $secondScore)  = split('\s+', $line); 
		    } #End peptide processing
		    
		    last; #we have finished the peptide.  Theoretically if the protein has more than one cleavage site we could keep going but will just process this sequence again instead
		}
	    }
	    
	    next if $skipProcessing; #proceed to next peptide
	    
	    #add results to $peptideInfo
	    if ($disorderString eq ""){   #did not find expected peptide start position in disopred result file
		$self->writeError($modbaseSeqId, $peptideStartPosition, "N/A", "no_disopred_peptide_sequence", "processDisopredResults", 
				  "Found Disopred results file but did not find peptide sequence start position $peptideStartPosition in file", 1, 0);  #tested
		$errorCount++;
	    }
	    else {
		$peptideCount++;
		$peptideInfo->{disopredString} = $disorderString;
		my $firstScoreString = join (',', @firstScoreArray);
		$peptideInfo->{disopredFirstScores} = $firstScoreString;
	    }
	}
    }

    $self->writeLog("Done running Disopred. Processed $sequenceCount sequences containing $peptideCount peptides. Encountered $errorCount errors during this run");
    
    my $endTime = time();
    $self->writeDuration($startTime, $endTime, "runDisopred");
}




############################################################################################################################################################      
# runDisopredCommand
# Run the disopred executable to predict disorder for the given modbase sequence id. After running disopred, copy results to the global data directory.
# Two files are generated: (1) $modbaseSequenceId.diso and (2) $modbaseSequenceId.horiz_d. The first has the actual calls and the second is for reference. 
#
# PARAM  $modbaseSequenceId: sequence id on which to run disopred
# PARAM  $peptideStartPostion: peptide that is being processed when this call is made (only needed for error handling)
# PARAM  $disopredOutputDir: global data directory that stores all of the disopred results files (partitioned into three letter prefixe directories for sequences).
# RETURN -1 if an error occurs; any output from running Disopred otherwise
#############################################################################################################################################################  
sub runDisopredCommand{
    my ($self, $modbaseSequenceId, $peptideStartPosition, $disopredOutputDir) = @_;
    

    #create input fasta file (one sequence only)
    my $runDisopredCmd = $self->getParam("run_disopred_cmd");

    my $cwd = getcwd();
    
    my $tempDir = File::Temp->tempdir("XXXX", CLEANUP =>1);
    chdir($tempDir);

    my $tempFastaFile = "$modbaseSequenceId.fasta";
    my $fh = FileHandle->new(">" . $tempFastaFile);
    unless ($fh =~ /FileHandle/){
	$self->writeError("$modbaseSequenceId", $peptideStartPosition, "N/A", "file_missing", "runDisopredCommand", "could not open temporary fasta file $tempFastaFile for writing: $!", 1, 0);
	chdir($cwd);
	return -1;   #tested, but not present in regression test (couldn't figure out a way to break it. Make fake directory in $tempFastaFile and run regression test if further testing required).
    }
    my $sequence = $self->{AllSequences}->{$modbaseSequenceId}->{fullSequence};
    print $fh ">$modbaseSequenceId|$peptideStartPosition\n$sequence";

    #TODO -- see if we need to gather $shortSeqId or if we can just use $modbaseSequenceId

    my $shortSeqId = $tempFastaFile;
    if ($tempFastaFile =~ /(.*)\.fasta/){
	$shortSeqId = $1;
    }

    my $cmd = "$runDisopredCmd $tempFastaFile";
    $self->writeLog("running disopred command $cmd");
    my $disopredOutput = $self->runSimpleCommand($cmd, 0);

    #copy results to disopred global results directory (will read them from there the next time it's needed for this sequence)
    my $mkdirCmd = "mkdir -p $disopredOutputDir";
    my $output = $self->runSimpleCommand($mkdirCmd, 1);

    my $horizCpCmd = "cp $shortSeqId.horiz_d $disopredOutputDir/$modbaseSequenceId.horiz_d";
    $output = $self->runSimpleCommand($horizCpCmd, 1);

    my $disoCpCmd = "cp $shortSeqId.diso $disopredOutputDir/$modbaseSequenceId.diso";
    $output = $self->runSimpleCommand($disoCpCmd, 1);

    chdir($cwd);
    unlink $tempDir;
    return $disopredOutput;

}


#############################################################################################################################################################  
# runSimpleCommand
# Runs system command and logs output if requested.
# 
# PARAM  $cmd: the command to run
# PARAM  $writeOutput: whether to log output if its present
# RETURN $output: command stderr / stdout
#############################################################################################################################################################  
sub runSimpleCommand{
    my ($self, $cmd, $writeOutput) = @_;
    my $output = `$cmd 2>&1`;
    if ($output && $writeOutput){
	$self->writeLog("ran command $cmd; exited with non-zero output $output");
    }
    return $output;
}

#############################################################################################################################################################  
# makePsipredStructureMap
# Simple method that maps psipred calls to normal structure types
# RETURN $map: dictionary where keys are psipred calls and values are the structure types
#############################################################################################################################################################  
sub makePsipredStructureMap{
    my ($self) = @_;
    my $map;
    $map->{"C"} = "L";
    $map->{"E"} = "B";
    $map->{"H"} = "A";
    return $map;

}


############################################################################################################################################################      
# runPsipredCommand
# Run the psipred executable to predict secondary structure for the given modbase sequence id. After running psipred, copy results to the global data directory.
# Two files are generated: (1) $modbaseSequenceId.ss2 and (2) $modbaseSequenceId.horiz. The first has the actual calls and the second is for reference. 
#
# PARAM  $modbaseSequenceId: sequence id on which to run psipred
# PARAM  $peptideStartPostion: peptide that is being processed when this call is made (only needed for error handling)
# PARAM  $psipredOutputDir: global data directory that stores all of the psipred results files (partitioned into three letter prefix directories for sequences).
# RETURN -1 if an error occurs; any output from running PsiPred otherwise
#############################################################################################################################################################  
sub runPsipredCommand{
    my ($self, $modbaseSequenceId, $peptideStartPosition, $psipredOutputDir) = @_;
    
    #create input fasta file (one sequence only)
    my $runPsipredCmd = $self->getParam("run_psipred_cmd");

    my $cwd = getcwd();
    
    my $tempDir = File::Temp->tempdir("XXXX", CLEANUP =>1);
    chdir($tempDir);


    my $tempFastaFile = "$modbaseSequenceId.fasta";

    my $fh = FileHandle->new(">" . $tempFastaFile);

    unless ($fh =~ /FileHandle/){
	$self->writeError("$modbaseSequenceId", $peptideStartPosition, "N/A", "file_missing", "runPsipredCommand", "could not open temporary fasta file $tempFastaFile for writing: $!", 1, 0);
	chdir($cwd);
	return -1;   #tested, but not present in regression test (couldn't figure out a way to break it. Make fake directory in $tempFastaFile and run regression test if further testing required).
    }
    my $sequence = $self->{AllSequences}->{$modbaseSequenceId}->{fullSequence};
    print $fh ">$modbaseSequenceId|$peptideStartPosition\n$sequence";

    #TODO -- see if we need to gather $shortSeqId or if we can just use $modbaseSequenceId

    my $shortSeqId = $tempFastaFile;
    if ($tempFastaFile =~ /(.*)\.fasta/){
	$shortSeqId = $1;
    }

    my $cmd = "$runPsipredCmd $tempFastaFile";
    $self->writeLog("running psipred command $cmd");
    my $psipredOutput = $self->runSimpleCommand($cmd, 0);

    #copy results to psipred global results directory (will read them from there the next time it's needed for this sequence)
    my $mkdirCmd = "mkdir -p $psipredOutputDir";
    my $output = $self->runSimpleCommand($mkdirCmd, 1);
    
    my $horizCpCmd = "cp $shortSeqId.horiz $psipredOutputDir/$modbaseSequenceId.horiz";
    $output = $self->runSimpleCommand($horizCpCmd, 1);

    my $psipredCpCmd = "cp $shortSeqId.ss2 $psipredOutputDir/$modbaseSequenceId.ss2";
    $output = $self->runSimpleCommand($psipredCpCmd, 1);

    chdir($cwd);
    unlink $tempDir;

    return $psipredOutput;

}


#############################################################################################################################################################  
# getProteinNames
# Read file containing uniprot accessions and their names, and add this annotation to the protein sequences currently being processed. File is already generated
# (currently from full uniprot annotation file). File only contains human names; consider expanding to include other species.
#
# RETURN NULL
#############################################################################################################################################################  
sub getProteinNames{
    my ($self) = @_;
    
    my $proteinNameFile = $self->getParam("protein_name_file");

    my $proteinNameFh = FileHandle->new("<" . $proteinNameFile);
    unless ($proteinNameFh =~ /FileHandle/){

	$self->writeError("global", "global", "global", "file_missing", "getProteinNames", "Could not open protein name file $proteinNameFile: $!", 1, 1);
	return 0;
    }

    my $proteinNames;

    while (<$proteinNameFh>){
	chomp;
	my $line = $_;
	my ($accession, $proteinName) = split('\t', $line);

	$proteinNames->{$accession} = $proteinName;
    }

    my $allPeptides = $self->{AllSequences};
    foreach my $modbaseSeqId (keys %$allPeptides){
	my $uniprotAccession = $allPeptides->{$modbaseSeqId}->{uniprotAccession};
	my $proteinName = $proteinNames->{$uniprotAccession};

	$self->{AllSequences}->{$modbaseSeqId}->{proteinName} = $proteinName;
    }
}



#############################################################################################################################################################  
# makeColumnHeaderString
# Reads $self->{ColumnInfo} and gets all column names that will be output in the results file. Creates a tab-delimited string with the names. Columns are ordered
# according to specification in $self->getParam("column_info_file") which was loaded previously.
#
# RETURN $columnHeaderString: tab-delimited string with all names of columns to be output.
#############################################################################################################################################################  
sub makeColumnHeaderString{

    my ($self) = @_;

    my $columnHeaderString = "";
    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");

    my $columns = $self->{ColumnInfo};
    unless ($columns){
	$self->writeOutputError("Columns were not loaded correctly (did not find anything in self->{ColumnInfo})"); #tested
    }
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   
	
	#check PeptidePipeline.pm outputs this column
 	if ($modes =~ /$keywordFeatureInternal/){                          
	    my $displayName = $columns->{$columnShortName}->{displayName};
	
	    #check we have run the method that would generate this column
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


#############################################################################################################################################################  
# printAllPeptides
# Simple method that prints all peptides and their features to STDERR. Not called by anything else in this module,  but useful as a convenience function.
#
# RETURN Null
#############################################################################################################################################################  
sub printAllPeptides{
    my ($self) = @_;
    return 0;

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


#############################################################################################################################################################  
# makeOutputLine
# Given all feature names and values for a particular peptide, create a tab-delimited string with the values. 
# Uses $self->{ColumnInfo} to coordinate order and which values to write.
# 
# RETURN Null
#############################################################################################################################################################  
sub makeOutputLine{
    my ($self, $outputInfo) = @_;
    my $columns = $self->{ColumnInfo};
    my $outputLine = "";
    my $keywordFeatureInternal = $self->getParam("keyword_feature_internal");
    foreach my $columnShortName (sort ({$columns->{$a}->{displayOrder} <=> $columns->{$b}->{displayOrder}} keys %$columns)){
	my $modes = $columns->{$columnShortName}->{modes};   
	
	#check PeptidePipeline.pm outputs this column
	if ($modes =~ /$keywordFeatureInternal/){

	    #check we have run the method that would generate this column
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


#############################################################################################################################################################  
# getResidueSolventAcc
# Get total possible exposed surface area by residue.
# Taken from Rose et al, Science, 1985, which determined these values using atomic radii calculated in Lee & Richards, 1971
#
# RETURN $residues: hash where keys are one-letter residue codes and values are possible surface area for that residue, in squared angstroms
#############################################################################################################################################################  
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

#############################################################################################################################################################  
# getSolventExposureFraction
# Calculate the fraction of the given residue's surface area that is exposed to solvent, and also our call for whether this is solvent-exposed or not.
# Exposure is anything that is > .033 fraction exposed.
#
# PARAM $accessibility: the observed exposed surface area of the residue, according to DSSP
# PARAM $residueOneLetter: one-letter code for the residue in question
# PARAM $residueSolventAcc: hash where the keys are the one letter residue codes and the values are the total possible exposed surface area for that residue
#
# RETURN Two element array; first entry is the fraction exposed and the second is the call ("A" for accessible and "N" for non-accessible)
#############################################################################################################################################################  
sub getSolventExposureFraction{
    my ($self, $accessibility, $residueOneLetter, $residueSolventAcc) = @_;

    $accessibility *= 1.0;
    my $totalSa = $residueSolventAcc->{$residueOneLetter};
    return (0, 0) unless $totalSa;
   
    my $fractionSa = $accessibility / $totalSa;
   
    my $accCall;
    if ($fractionSa > 0.33){
	$accCall .= "A";
    }
    else {
	$accCall .= "N";
    }

    return ($fractionSa, $accCall);
}


#############################################################################################################################################################  
# writeUserResults
# Writes all peptide features to the pipeline results file. Features are all stored in hash at $self->{AllSequences}->{$modbaseSeqId}->{$peptideStartPosition}
# One tab-delimited feature line is written for each peptide.  Some features are always written and some depend on whether the method was run to generate them. Uses $self->{ColumnInfo}
# to coordinate this (along with order in which features are written).
#
# RETURN Null
#############################################################################################################################################################  
sub writeUserResults{

    my ($self) = @_;

    my $startTime = time();
    $self->writeLog("Writing results");

    return unless $self->runMethod("writeUserResults");
    
    my $resultsFh = $self->{ResultFh};
    my $sequenceCount = 0;
    my $peptideCount = 0;
    my $errorCount = 0;
    
    #write column header line
    eval {
	my $columnHeaderString = $self->makeColumnHeaderString();
	
	print $resultsFh $columnHeaderString . "\n";
	
	my $residueSolventAcc = $self->getResidueSolventAcc();
	
	my $allSequences = $self->{AllSequences};
	
	return 0 if ($self->processGlobalErrors());
	
	foreach my $modbaseSeqId (keys %$allSequences){
	    $sequenceCount++;
	    my $peptides = $allSequences->{$modbaseSeqId}->{peptides};
	    
	    my $uniprotAccession = $allSequences->{$modbaseSeqId}->{uniprotAccession};
	    my $proteinName = $allSequences->{$modbaseSeqId}->{proteinName};

	    my $totalModelCount = $allSequences->{$modbaseSeqId}->{totalModelCount};
	    
	    foreach my $peptideStartPosition (keys %$peptides){
		
		$peptideCount++;
		my $peptideInfo = $peptides->{$peptideStartPosition};
		my $bestModelId = $peptideInfo->{bestModelId};
		
		#outputInfo is the final hash containing output values; keys are feature names and values are their values.
		my $outputInfo;
		
		#write general peptide information    
		my @attributes = ("modelsContainingPeptideCount", "run",  "modelTargetStart", "modelTargetEnd", "peptideEndPosition", "peptideSequence",  "modelScore", "templateSequenceIdentity", "coverage", "loopLength",  "classification",  "disopredString", "disopredFirstScores",  "psipredString", "psipredFirstScores",  "fullAlignmentFilePath");
		
		foreach my $attribute (@attributes){
		    my $value = $peptideInfo->{$attribute};
		    $outputInfo = $self->addOutputValue($outputInfo, $value, $attribute);
		}
		$self->addOutputValue($outputInfo, $peptideStartPosition, "peptideStartPosition");  #get this from $peptide key instead of feature (helps with errors)
		
		my $url = "";
		if ($bestModelId){
		    $url = "http://salilab.org/modbase/search?modelID=" . $bestModelId . "&displaymode=moddetail";
		}
		
		#write data that is common to all peptides for this $modbaseSeqId
		$outputInfo = $self->addOutputValue($outputInfo, $modbaseSeqId, "sequenceId");
		$outputInfo = $self->addOutputValue($outputInfo, $uniprotAccession, "uniprotAccession");
		$outputInfo = $self->addOutputValue($outputInfo, $proteinName, "proteinName");
		$outputInfo = $self->addOutputValue($outputInfo, $totalModelCount, "totalModelCount");
		$outputInfo = $self->addOutputValue($outputInfo, $bestModelId, "modelId");
		$outputInfo = $self->addOutputValue($outputInfo, $url, "url");
		
		#if there was at least one model containing the peptide, write secondary structure, solvent accessibility, and data regarding model details
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
		
		#Any errors that occurred
		my $errors = $peptideInfo->{errors};
		my $errorString = "";
		if ($errors){
		    $errorCount++;
		    $errorString = join (',', @$errors);
		}
		else {
		    $errorString = $self->getParam("keyword_no_cluster_errors");
		}
		
		$outputInfo = $self->addOutputValue($outputInfo, $errorString, "errors");
		
		#write final line for this peptide
		my $outputLine = $self->makeOutputLine($outputInfo);
		print $resultsFh $outputLine . "\n";
	    }
	}
	my $noPeptidesParsedSeqs = $self->getNoPeptidesParsedSeqs();
	foreach my $seqId (keys %$noPeptidesParsedSeqs){
	    my $noPeptidesParsedKeyword = $self->getParam("keyword_no_peptides_parsed");
	    my $uniprotAccession = $allSequences->{$seqId}->{uniprotAccession};
	    my $outputInfo;
	    $outputInfo = $self->addOutputValue($outputInfo, $seqId, "sequenceId");
	    $outputInfo = $self->addOutputValue($outputInfo, $uniprotAccession, "uniprotAccession");
	    $outputInfo = $self->addOutputValue($outputInfo, $noPeptidesParsedKeyword, "errors");
	    my $outputLine = $self->makeOutputLine($outputInfo);
	    print $resultsFh $outputLine . "\n";
	}
	    
    };
    unless ($@){
	
	$self->writeLog("Done writing results. Wrote lines for $sequenceCount sequences containing $peptideCount peptides. $errorCount peptides had at least one error");
	
	my $endTime = time();
	$self->writeDuration($startTime, $endTime, "writeUserResults");
    }
}


#############################################################################################################################################################  
# addOutputValue
# Simple method that adds a feature name / value pair to a hash, after checking that the feature is a valid feature name.
#
# PARAM  $outputInfo: Hash where keys are feature names and values are their values
# PARAM  $value: feature value
# PARAM  $internalName: feature name
# RETURN $outputInfo: same hash as input after it has been modified.
#############################################################################################################################################################  
sub addOutputValue{
    my ($self, $outputInfo, $value, $internalName) = @_;
    
    my $columns = $self->{ColumnInfo};
    my $columnInfo = $columns->{$internalName};
    
    #check valid column name. 
    unless ($columnInfo){
	my $keyString = join(", ", keys %$columns);
	my $msg = "trying to add value to column info with internal id $internalName which is not a defined internal id.  List of these IDs: $keyString\n";  
	$self->writeOutputError($msg); #tested
    }
    $outputInfo->{$internalName} = $value;

    return $outputInfo;
}



sub writeOutputError{
    my ($self, $msg) = @_;
    $msg = "ERROR: $msg";
    $self->writeLog($msg);
    $self->writeResultFhOutputError();
    
    $self->finalize();
    die "$msg"; # not trapping message since it is written to log here
    
}


#############################################################################################################################################################  
# processGlobalErrors
# If pipeline had massive global errors (such as not finding an internal file that is supposed to exist), write this to the Results file along with the error code.
# This is intended to bypass all other output.
#
# RETURN 1 if there was a global error; 0 otherwise
#############################################################################################################################################################  
sub processGlobalErrors{
    my ($self) = @_;
    
    if ($self->pipelineHasErrors()){  #only set if there is a 'global' error
	my $sequences = $self->{AllSequences};
	my $errors = $sequences->{"global"}->{peptides}->{"global"}->{errors};  #hackish way of representing errors
	my $errorString;
	if ($errors){
	    $errorString = join (',', @$errors);
	}
	else {
	    $self->writeOutputError("found global error but no error code written");  
	    #tested but not in the regression test; too complicated to make an error case for this, and also shouldn't happen unless $allSequences has a key error, 
	    #which isn't being tested for anywhere else, for the most part
	}
	my $outputInfo;
	$outputInfo = $self->addOutputValue($outputInfo, $errorString, "errors");
	$outputInfo = $self->addOutputValue($outputInfo, "global", "sequenceId");	
	$outputInfo = $self->addOutputValue($outputInfo, "global", "peptideStartPosition");
		
	my $outputLine = $self->makeOutputLine($outputInfo);
	
	my $resultFh = $self->{ResultFh};
	print $resultFh $outputLine . "\n";
	$self->writeLog("Done writing results. Wrote lines for 1 sequences containing 1. 1 peptides had at least one error");
	return 1;
    }
    else {
	return 0;
    }
}



#############################################################################################################################################################  
# parseAlignmentFile
# Read alignment file, determine which entry is the target and which is the template (varies according to modpipe run), and return info on both.
# 
# PARAM  $fullAlignmentFileName: full path and name of alignment file
# PARAM  $runName: name of modpipe run that generated this alignment
#
# RETURN Four element array: (1) PIR header line of target (2) target protein sequence (3) PIR header line of template (4) template protein sequence
#                                PIR format description: http://salilab.org/modeller/manual/node454.html 
#############################################################################################################################################################  
sub parseAlignmentFile{

    my ($self, $fullAlignmentFileName, $modpipeRunName, $modbaseSeqId, $peptideStartPosition, $modelId) = @_;
    
    my $runInfo = $self->getRunInfo();
    my $runStyle = $runInfo->{$modpipeRunName}->{style};
    
    my $alignmentFh = FileHandle->new("<" . $fullAlignmentFileName);
    unless ($alignmentFh =~ /FileHandle/){
	$self->writeError($modbaseSeqId, $peptideStartPosition, $modelId, "file_missing", "parseAlignmentFile", "Could not open alignment file $fullAlignmentFileName: $!", 1, 0);
	return 0; #Untested (after $runInfo parameterized)
    }
	
    
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
    #Ali file has two proteins (template and target). New-style modpipe runs have target as first protein; old style has target as second protein
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
	$self->writeError($modbaseSeqId, $peptideStartPosition, $modelId, "internal_error",  "parseAlignmentFile",   #TODO - consider making modpipe_run_info_error after adding runinfo parameters
			  "Did not get expected modpipe run style for run $modpipeRunName: expected 'old' or 'new', got $runStyle", 1, 0); #Untested (after $runInfo parameterized)
	return 0;
    }

    return ($targetInfoLine, $targetSequence, $templateInfoLine, $templateSequence);

}


#############################################################################################################################################################  
# getModelOffset
# Some runs in Modpipe output a model file that begins counting residues with 1 instead of the index of the first modeled residue in the protein. To account for 
# this when searching for our peptide in the model, we need to subtract the index of the first modeled residue from the peptide's position. This method returns
# that first modeled residue, known as the offset, after checking whether it's necessary by examining which modpipe run generated the model.
#
# PARAM  $peptideInfo: normal $peptideInfo hash (format described in parsePeptidesFromSequences)
# RETURN $offset: value to subtract from peptide range; may be 0 depending on the modpipe run.
#############################################################################################################################################################  
sub getModelOffset{
 
    my ($self, $peptideInfo) = @_;

    #get run info, validate
    my $run = $self->getParam("modpipe_run_name");   
    if (($run ne "human") && ($run ne "human_4-2007") && ($run ne "snp-human2")  && ($run ne "snp-human3")  && ($run ne "snp-human4") && ($run ne "trembl2004") && ($run ne "human_2008" )
	&& ($run ne "pcss_modpipe_runs")){
	return "Did not get valid modpipe run name (one of 'human', 'human_4-2007', 'snp-human2', 'snp-human3', 'snp-human4', 'termbl2004', 'human_2008', 'pcss_modpipe_runs': instead got $run";
    }

    #use run to check if there should be offset
    my $runInfo = $self->getRunInfo();
    my $hasOffset = $runInfo->{$run}->{offset};
    my $offset;

    if ($hasOffset){
	$offset = $peptideInfo->{modelTargetStart} - 1;  #offset is equal to the first position in the protein that is in the model (-1, as model numbering is base-1)
    }
    else {
	$offset = 0;
    }
    return $offset;
}



#############################################################################################################################################################  
# getModelPeptideRange
# Given a peptide range, return position of the peptide as it appears in the model (may change due to certain Modpipe runs renumbering residues in the models;
# see getModelOffset() for details.
# 
# PARAM  $peptideInfo:  normal $peptideInfo hash (format described in parsePeptidesFromSequences)  
# RETURN Two element array; entries are the start and end residue positions of the peptide in the model file.
#############################################################################################################################################################  
sub getModelPeptideRange{

    my ($self, $peptideInfo) = @_;

    my $offsets;

    my $offset = $self->getModelOffset($peptideInfo);  #check if peptide range should be offset if we are searching for it in the model

    if ($offset =~ /\D/){  #error message
	return $offset;
    }

    #no errors; adjust peptide ranges and return
    my $peptideStartPosition = $peptideInfo->{peptideStartPosition};
    my $peptideEndPosition = $peptideInfo->{peptideEndPosition};
    
    $peptideStartPosition -= $offset;
    $peptideEndPosition -= $offset;

    return ($peptideStartPosition, $peptideEndPosition);
}


#############################################################################################################################################################  
# checkInput
# After sequences have been loaded, call this method to make sure that at least one peptide is present and will be processed. This might not be the case
# in Application Scan mode where the supplied sequences have been parsed but no peptides matching the user-supplied rules have been discovered. This method
# will indicate such and the rest of the pipeline does not have to proceed.
#
# RETURN 1 if there is at least one peptide, 0 otherwise
#############################################################################################################################################################  
sub checkInput{
    my ($self) = @_;
    my $allSequences = $self->{AllSequences};

    foreach my $seqId (keys %$allSequences){
	my $peptides = $allSequences->{$seqId}->{peptides};
	my $peptideCount = scalar(keys %$peptides);
	if ($peptideCount > 0){
	    return 1;
	}
    }
    return 0;
}


#############################################################################################################################################################  
# writeNoInput
# The equivalent of writeUserResults() except if there were no peptides to be processed. Writes the column header and then indicates that no peptides were parsed.
# 
# RETURN Null
#############################################################################################################################################################  
sub writeNoInput{
    my ($self) = @_;

    my $resultsFh = $self->{ResultFh};

    my $columnHeaderString = "";
    
    eval {

	$columnHeaderString = $self->makeColumnHeaderString();  
	
    };

    unless ($@){

	print $resultsFh $columnHeaderString . "\n";
	my $noPeptidesParsedKeyword = $self->getParam("keyword_no_peptides_parsed");
	print $resultsFh $noPeptidesParsedKeyword . "\n";
    }
}

#############################################################################################################################################################  
# loadErrorCodes
# Initializes error code hash. These is an attempt at a CV of possible errors, and a longer description for each (values of the hash). Methods that have an error
# will write which one it was by referring to one of these codes. Note that "internal_error" MUST be kept as one of these keys; otherwise calling $self->writeError()
# and passing an invalid error code will result in an infinite loop of never being able to find the appropriate error code (this is checked here).
#
# RETURN Null
#############################################################################################################################################################  
sub loadErrorCodes{

    my ($self) = @_;
    $self->{ErrorCodes}->{"no_dssp_structure_info"} = "Sequences for which there were models but DSSP secondary structure output file could not be read at the cleavage sequence position.";
    $self->{ErrorCodes}->{"no_model_in_expected_directory"} = "Sequences for which there was a best model but it was not found in expected source modbase directory.";
    $self->{ErrorCodes}->{"dssp_format_error"} = "DSSP files with unexpected format for a given line";

    $self->{ErrorCodes}->{"no_disopred_result_file"} = "Sequences for which no Disopred results file was found.";
    $self->{ErrorCodes}->{"disopred_sequence_mismatch"} = "Sequences for which the residues at the cleavage sequence positions in the Disopred results file did not match up with the annotated cleavage sequence.";
    $self->{ErrorCodes}->{"no_disopred_peptide_sequence"} = "Sequences for which there was a DISOPRED results file but the start position of the cleavage sequence was not found in the file.";
    $self->{ErrorCodes}->{"disopred_format_error"} = "Disopred files with unexpected format for a given line";

    $self->{ErrorCodes}->{"no_prof_result_file"} = "Sequences for which no Prof results file was found.";
    $self->{ErrorCodes}->{"prof_sequence_mismatch"} = "Sequences for which the residues at the cleavage sequence positions in the Prof results file did not match up with the annotated cleavage sequence.";

    $self->{ErrorCodes}->{"no_psipred_result_file"} = "Sequences for which no Psipred results file was found.";
    $self->{ErrorCodes}->{"psipred_sequence_mismatch"} = "Sequences for which the residues at the cleavage sequence positions in the Psipred results file did not match up with the annotated cleavage sequence.";
    $self->{ErrorCodes}->{"no_psipred_peptide_sequence"} = "Sequences for which there was a PSIPRED results file but the start position of the cleavage sequence was not found in the file.";
    $self->{ErrorCodes}->{"psipred_format_error"} = "Psipred files with unexpected format for a given line";

    $self->{ErrorCodes}->{"unexpected_pdb_code_in_alignment"} = "Alignment files that did not have a four letter pdb code in the expected column in the pir line for the template.";
    $self->{ErrorCodes}->{"incorrect_template_position_in_alignment"} = "Alignment files that appear to have the template in the wrong position (first or second relative to the target) as expected according to the modpipe style";

    $self->{ErrorCodes}->{"mismatch_cleavage_sequence_in_alignment_file"} = "Alignment files that did not have the annotated cleavage sequence at the expected position for the target sequence.";

    $self->{ErrorCodes}->{"unexpected_lysate_call"} = "Entries in lysate file that did not have 'A' or 'P' in expected column";

    $self->{ErrorCodes}->{"file_missing"} = "File that necessary for global peptide operation (eg sequence input, rules file, etc.)";

    $self->{ErrorCodes}->{"no_defined_peptides"} = "Modbase sequence did not have any peptides listed in fasta header, or only mismatches were listed";

    $self->{ErrorCodes}->{"invalid_model"} = "Models that are invalid models (errors in sequence ranges, etc.)";
    $self->{ErrorCodes}->{"regex_error"} = "Files which were parsed expecting (but not getting) a regular expression";
    $self->{ErrorCodes}->{"column_error"} = "Columns were not loaded correctly";
    $self->{ErrorCodes}->{"internal_error"} = "General internal error not dependent on any input issues";    
    $self->{ErrorCodes}->{"modpipe_run_info_error"} = "Errors with modpipe run info data structure";

    unless ($self->{ErrorCodes}->{"internal_error"}){
	$self->writeLog("internal_error MUST be a defined error code to avoid infinte recursive error loops. Please set 'internal_error' in loadErrorCodes().");
	exit(1);  
    }

}


#############################################################################################################################################################  
# checkResidueMatchesDssp
# Quick check to make sure the residues in the given peptide match up with their corresponding residues in the DSSP file.
# 
# PARAM  $dsspResidueOneLetter:  the residue as it appears in the dssp file
# PARAM  $position in peptide: the position of the residue in the peptide (1-based)
# PARAM  $peptideSequence: full peptide sequence
#
# RETURN: 1 if there is a match; an error message otherwise
#############################################################################################################################################################  
sub checkResidueMatchesDssp{

    my ($self, $dsspResidueOneLetter, $positionInPeptide, $peptideSequence) = @_;

    my @residues = split('', $peptideSequence);
    my $residueInPeptide = $residues[$positionInPeptide - 1];  #TODO - after testing implemented, change $positionInPeptide to be zero-based

    if ($dsspResidueOneLetter ne $residueInPeptide){
	my $msg = "Error: peptide has residue $residueInPeptide at position $positionInPeptide but in dssp file the residue is $dsspResidueOneLetter";
	return $msg;
    }
    else {
	return 1;
    }
}

#############################################################################################################################################################  
# writeError
# If an error occurs in the pipeline, call this method to handle it. Takes an error code which must be validated as an approved error code (so the server backend
# can handle it if necessary). Errors can either be global (generally catastrophic, such as when an expected file is missing) or "feature" which means that 
# one feature of the modbase sequence and peptide id combination is missing, but other features can be processed. 
#
# If the error is a global error, an entry "global" is made in $self->{AllSequences} with one peptide in its peptides hash (also named "global"). If it is a feature
# error, then the error is assigned to the sequence and peptide in question. Either way, the error goes into an array for the 'sequence/peptide' combo (keyed on "errors").
# Global errors will set {SkipRemainingPipeline} which other top-level pipeline methods check to see if they should continue. The error is logged and also written to
# the results file "errors" column for the sequence / peptide combination (if it is a global error, then this should be the only line in the results file, with only
# the sequence, peptide, and errors columns filled).
#
# PARAM  $modbaseSequenceId, $peptideId: modbase sequence id and peptide start position of the peptide that threw the error (or "global" if it is a global error).
# PARAM  $modelId: if the error was the result of trying to process a model, this is also passed (usually "N/A" if not).
# PARAM  $errorCode: error code for the error from a controlled vocab; specified in $self->loadErrorCodes().
# PARAM  $methodName: name of the method that threw the error, for reference (not checked).
# PARAM  $errorMessage: well-crafted and specific error message
# PARAM  $writeLog: whether to log the error mesage or not
# PARAM  $skipRemainingPipeline: whether to skip the rest of the pipeline or not (always the case with global errors, although not checked).
# RETURN NULL
#############################################################################################################################################################  
sub writeError{

    my ($self, $modbaseSeqId, $peptideId, $modelId, $errorCode, $methodName, $errorMessage, $writeLog, $skipRemainingPipeline) = @_;

    my $errorDescription = $self->{ErrorCodes}->{$errorCode};
    unless ($errorDescription){
	my $errorCodes = $self->{ErrorCodes};
	my $errorList = join (", ", keys %$errorCodes);
	my $msg = "Attempted to write an error code ($errorCode) that is not a defined error code. List of allowed codes: $errorList";
	
	$self->writeError("global", "global", "global", "internal_error", "writeError", "$msg", 1, 1);   
	#note this results in an infinite loop if "internal_error" isn't set. There is a check to ensure "internal_error" is set upon initial call to loadErrorCodes()
	
	return 0;  
	#Note that since writeError is called from so many places, there is a good chance that processing will continue in that method. This should be fine although I haven't gone through
	#and tested to make sure that we get expected output if someone writes a bad error code. After that method ends, pipeline will skip to writeUserResults() due to having a global error
	#set here and will proceed with normal global error processing 
    }
    push (@{$self->{AllSequences}->{$modbaseSeqId}->{peptides}->{$peptideId}->{errors}}, $errorCode);

    if ($skipRemainingPipeline){
	$self->{SkipRemainingPipeline} = 1;
    }

    if ($writeLog){
	$self->writeLog("ERROR:  $methodName modbase seq id: $modbaseSeqId peptide start postion: $peptideId model id: $modelId Skip remaining pipeline: $skipRemainingPipeline Error: $errorMessage");
    }
}



#############################################################################################################################################################  
# getRunInfo
# Get {RunInfo} hash that contains all details for modpipe run that created the models to be examined for structural features.
#
# RETURN Run Info hash (see loadRunInfo() for specifics
#############################################################################################################################################################  
sub getRunInfo{
    my ($self) = @_;
    return $self->{RunInfo};
}

#############################################################################################################################################################  
# loadRunInfo
# Load necessary information for different modpipe runs that allows us to find alignment and model paths, etc.
# Each entry contains the following:
# {path} = full path of alignment and model top level directories
# {offset} = 1 if the residue position in the model file has to be adjusted to match the position of the peptide that's being processed  (see $getModelOffset() for details)
# {style} = modpipe run style; the most important thing is specifies is if the template or target is the first entry in an alignment file.
# 
# This method could also be loaded from a file so that it's not hardcoded (on the todo list).
#
# RETURN NULL
#############################################################################################################################################################  
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

    $runInfo->{"pcss_modpipe_runs"}->{path} = "/netapp/sali/peptide/data/landing/alignments/";
    $runInfo->{"pcss_modpipe_runs"}->{offset} = 0;
    $runInfo->{"pcss_modpipe_runs"}->{style} = "new";    



    $runInfo->{"trinidad"}->{path} = "/park2/modbase/projects/trinidad/data";
    $runInfo->{"trinidad"}->{offest} = 0;
    $runInfo->{"trinidad"}->{style} = "new";

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


#############################################################################################################################################################  
# Return 1 if pipeline has global or fatal errors; intended to be called by top level methods before they run.
#############################################################################################################################################################  
sub pipelineHasErrors{
    my ($self) = @_;
    if ($self->{SkipRemainingPipeline}){
	$self->writeLog("Skipping method as fatal errors were encountered earlier in the pipeline");
	return 1;
    }
    else {
	return 0;
    }
}


#############################################################################################################################################################  
# fulfillsMethodDependencies
# This, and other 'method dependency' functions are legacy when one method had to be run before another one (ie getBestModels() before parseDsspResults().
# Not used any more and could be deleted, although might be useful for extreme error checking.
#############################################################################################################################################################  
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

#############################################################################################################################################################  
# Adds a sequence to the list of those for which no peptides were parsed in Application Scan mode (stored as hash).
# 
# PARAM  $modbaseSeqId: Sequence ID to add
# RETURN NULL
#############################################################################################################################################################  
sub writeNoPeptidesParsedForSeq{
    my ($self, $modbaseSeqId) = @_;
    $self->{NoPeptidesParsedSeqs}->{$modbaseSeqId} = 1;
}

sub getNoPeptidesParsedSeqs{

    my ($self) = @_;
    my $noPeptidesParsedSeqs = $self->{NoPeptidesParsedSeqs};
    return $noPeptidesParsedSeqs;
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


return 1;


