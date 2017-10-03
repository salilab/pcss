include ../Makefile.include

.PHONY: install
D=${DATADIR}

FILES=${D}/caspaseBenchmarkScores_20100527 \
      ${D}/caspaseBenchmarkScores.txt \
      ${D}/caspaseCreateModel_svmModelFile ${D}/caspaseHomologMap.txt \
      ${D}/caspaseMissingOne ${D}/caspaseSvmModel_20100527 \
      ${D}/caspaseTrainingSetSvmModel.txt ${D}/columnInfo.txt \
      ${D}/globalPeptideParameters.txt ${D}/grbBenchmarkScores.txt \
      ${D}/grbCreateModel_svmModelFile ${D}/grbHomologMap.txt \
      ${D}/grb_homologs.txt ${D}/grbSvmModel_20100310.txt \
      ${D}/hivBenchmarkScores.txt ${D}/hivCreateModel_svmModelFile.txt \
      ${D}/human2008ModelTableTrinidad.txt ${D}/human2008ModelTable.txt \
      ${D}/human_2008_sql_results.rtf ${D}/leaveOneOutParams.txt \
      ${D}/outputFile.txt ${D}/params.txt \
      ${D}/peptidePipelineInput.txt ${D}/peptideRulesFile.txt \
      ${D}/testPeptideServerParameters.txt ${D}/trinidadRunInfo.txt \
      ${D}/uniprot100.txt ${D}/uniprot1.txt \
      ${D}/uniprotAccessionList.txt ${D}/uniprotDescriptions.txt \
      ${D}/uniprotO75791.txt

install: ${FILES}

${DATADIR}/%: %
	@if [ ! -d ${DATADIR} ]; then mkdir -p ${DATADIR}; fi
	install $< $@
