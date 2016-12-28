STRUCTURED_HOME = ..

include $(STRUCTURED_HOME)/Common/mk/Make.defs

lib_targets := $(shell ls src | \tr ' ' '\n' | $(SED) '/GNUmakefile/d' | sort | uniq)

#
# These "target-specific variable values" are the actions to perform on
# the sub directories
#   all         builds all the executables, including the libraries they depend
#               on
#   lib         builds all the libraries
#   test        builds all the test executables
#   run         runs all the test executables
#   clean       deletes files for this configuration
#
_all_actions = lib all test run clean
_all_subdir = src test
_action = lib
_subdir = src
lib       : _action = lib
lib       : _subdir = src
all       : _action = all
all       : _subdir = src test
test      : _action = test
test      : _subdir = test
run       : _action = run
run       : _subdir = test
clean     : _action = clean NODEPENDS=TRUE
clean     : _subdir = src test

.PHONY: $(_all_actions) $(lib_targets)

# This dependency ensures the action is done for every directory.
$(_all_actions) : $(lib_targets)

# Extra for cleaning
clean:
	-$(RM) *~ $(addprefix lib,$(addsuffix .$(_lib_ext),$(lib_targets)))

# This rule performs the action
$(lib_targets):
	@if [ "$(_action)x" = "libx" ]; then echo ""; echo "Checking $@ library..."; echo ""; fi
	@for x in $(_subdir); do \
           $(MAKE) --no-print-directory --directory=$$x/$@ $(_action); \
           done

printvars:
	@echo "lib_targets         : $(lib_targets)"
