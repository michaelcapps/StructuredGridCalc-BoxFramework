STRUCTURED_HOME = .

include $(STRUCTURED_HOME)/Common/mk/Make.defs

clean_targets := lib application/sandbox application/gpuSandbox \
        application/vlaloops application/laplacian \
        application/lapack application/cgnswrite application/latticeBoltzmann \
        application/wave

.PHONY: all doc clean $(clean_targets) dist

all:
	$(ECHO)echo "Enter a subdirectory to build an application, libraries, or tests.  Valid targets for this directory are 'doc' and 'clean'."

doc:
	$(mkdox)

dist: clean
	cd .. && tar --exclude=.hg* --exclude=doc/doxygen_sqlite3.db --exclude=doc/html* --exclude=doc/latex* -cvjf Structured.tar.bz2 Structured

clean: $(clean_targets)
	-$(RM) *~

# Rule for each clean target
$(clean_targets):
	$(MAKE) --no-print-directory --directory=$@ clean

printvars:
	@echo "clean_targets         : $(clean_targets)"
