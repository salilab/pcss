include ../Makefile.include

.PHONY: install

FILES=${DATA}/caspaseBenchmarkScores_20100527 \
      ${DATA}/caspaseBenchmarkScores.txt \
      ${DATA}/caspaseCreateModel_svmModelFile ${DATA}/caspaseHomologMap.txt \
      ${DATA}/caspaseMissingOne ${DATA}/caspaseSvmModel_20100527 \
      ${DATA}/caspaseTrainingSetSvmModel.txt ${DATA}/columnInfo.txt \
      ${DATA}/globalPeptideParameters.txt ${DATA}/grbBenchmarkScores.txt \
      ${DATA}/grbCreateModel_svmModelFile ${DATA}/grbHomologMap.txt \
      ${DATA}/grb_homologs.txt ${DATA}/grbSvmModel_20100310.txt \
      ${DATA}/hivBenchmarkScores.txt ${DATA}/hivCreateModel_svmModelFile.txt \
      ${DATA}/human2008ModelTableTrinidad.txt ${DATA}/human2008ModelTable.txt \
      ${DATA}/human_2008_sql_results.rtf ${DATA}/leaveOneOutParams.txt \
      ${DATA}/outputFile.txt ${DATA}/params.txt \
      ${DATA}/peptidePipelineInput.txt ${DATA}/peptideRulesFile.txt \
      ${DATA}/testPeptideServerParameters.txt ${DATA}/trinidadRunInfo.txt \
      ${DATA}/uniprot100.txt ${DATA}/uniprot1.txt \
      ${DATA}/uniprotAccessionList.txt ${DATA}/uniprotDescriptions.txt \
      ${DATA}/uniprotO75791.txt

install: ${FILES}

${DATADIR}/%: %
	@if [ ! -d ${DATADIR} ]; then mkdir -p ${DATADIR}; fi
	install $< $@
