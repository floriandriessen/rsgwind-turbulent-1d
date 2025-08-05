# Directory with extra modules for use in mod_usr.t
USER_SRCDIR := ../user-extra-modules

# Local user fortran files
FSRCS := mod_kind_parameter.f90
FOBJS := $(FSRCS:.f90=.o)
FMODS := $(FSRCS:.f90=.mod)

# Local AMRVAC files
TSRCS := mod_input_output.t
TOBJS := $(TSRCS:.t=.o)
TMODS := $(TSRCS:.t=.mod)

$(info Extra user source file(s) included: $(FSRCS) $(TSRCS))

# Build rule for local user object files
# Note amrvac's INC_DIRS because some user modules may depend on amrvac modules
%.o: $(USER_SRCDIR)/%.f90
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

# Build rule for local AMRVAC object files (.t -> .f conversion is auto done)
%.o: $(USER_SRCDIR)/%.f amrvac.h
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))

clean_user:
	echo 'Cleaning local user objects'
	$(RM) $(FOBJS) $(FMODS) $(TOBJS) $(TMODS)

# Rule for the targets
amrvac mod_usr.o: $(FOBJS) $(TOBJS)
