include ../Makefile.include

.PHONY: install

FILES=${BINDIR}/convertHuman2008results.pl ${BINDIR}/postprocessPeptideJob.pl \
      ${BINDIR}/postprocessTrainingJob.pl ${BINDIR}/preprocessPeptideJob.pl \
      ${BINDIR}/runModelPipeline.pl ${BINDIR}/runPeptidePipeline.pl

install: ${FILES}

${BINDIR}/%: %
	@if [ ! -d ${BINDIR} ]; then mkdir -p ${BINDIR}; fi
	install $< $@

