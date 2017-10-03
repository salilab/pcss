include Makefile.include

.PHONY: install

install:
	${MAKE} -C bin install
	${MAKE} -C lib install
	${MAKE} -C data install
