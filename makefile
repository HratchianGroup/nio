#
# NIO Makefile
#
# To use this makefile, one should set the FC macro to GNU or NVIDIA to
# identify which compiler type is being used. Note that the specific compiler
# version used when building NIO should match that used in building the module
# and library files provided by MQCPack.
#
# Typically, one's user environment includes $mqcinstall, which points to the
# location of the MQCPack library directories. Alternatively, the one can
# hardcode in the path to that directory in the line below that sets macro
# MQCDir.
#
# The macro LIBSALGEBRA lists how the LAPack and BLAS routines should be
# referenced by in the compiler invocation. They are set to the standard
# default values. If one wants to use alternate versions, one should modify the
# assignement to LIBSALGEBRA.
#
# No other modifications to the makefile should be needed.
#
#
# Set-up flags that may need changes based on user configuration...
#
  FC        = nvfortran
# FC        = gfortran
MQCDir      = $(mqcinstall)
LIBSALGEBRA = -llapack -lblas
USEOMP      = yes
#
# Set-up flags that should not need to be changed by the user...
#
ifeq ($(FC),gfortran)
  MQCMODS = $(MQCDir)/GNU/mod
  MQCLIB  = $(MQCDir)/GNU/lib
  ifeq ($(USEOMP),yes)
    OMPFLAGS = -fopenmp
  endif
  LIBS    =  $(LIBSALGEBRA) $(MQCLIB)/libmqc.a -lmkl_rt
  FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 $(OMPFLAGS)
else ifeq ($(FC),nvfortran)
  MQCMODS      = $(MQCDir)/NVidia/mod
  MQCLIB       = $(MQCDir)/NVidia/lib
  ifeq ($(USEOMP),yes)
    OMPFLAGS = -mp
  endif
  LIBS    = $(LIBSALGEBRA) -L$(MQCLIB)
  FCFLAGS = -Mallocatable=03 -r8 -i8 $(OMPFLAGS)
endif
#
# The 'all' rule.
#
all: nio.exe
#
# Generic rules for building module (*.mod) and object (*.o) files.
#
%.mod: %.F03
	$(FC) $(FCFLAGS) -I$(MQCMODS) -c $*.F03

%.o: %.F03
	$(FC) $(FCFLAGS) -I$(MQCMODS) -c $*.F03
#
# Generic rule for building general executable program (*.exe) from a standard
# f03 source (*.f03) file.
#
%.exe: %.F03 %_mod.F03 %_mod.mod $(MQCLIB)/libmqc.a
	$(FC) $(LIBS) $(Prof) -I$(MQCMODS) $(FCFLAGS) -o $*.exe $*.F03 $*_mod.o $(MQCLIB)/libmqc.a
#
# Clean rule for removing any object, module, or executable files in the working directory.
#
clean:
	rm -f *.exe *.mod *.o
