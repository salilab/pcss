include ../Makefile.include

.PHONY: install
PL=${PERLDIR}

FILES=${PL}/ApplicationPipeline.pm ${PL}/BenchmarkerPipeline.pm \
      ${PL}/Benchmarker.pm ${PL}/CleavageSiteModel.pm \
      ${PL}/CreationPipeline.pm ${PL}/FastaReader.pm \
      ${PL}/JackknifeBenchmarker.pm ${PL}/LeaveOneOutBenchmarker.pm \
      ${PL}/PeptidePipeline.pm ${PL}/Pipeline.pm ${PL}/SvmModel.pm

install: ${FILES}

${PL}/%: %
	@if [ ! -d ${PL} ]; then mkdir -p ${PL}; fi
	install $< $@

