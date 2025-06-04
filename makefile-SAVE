#
# This is a simple makefile for building spin-squared calculation code.
#
MQCDir       = $(mqcinstall)
MQCMODS      = $(MQCDir)/NVidia/mod
MQCLIB       = $(MQCDir)/NVidia/lib
LIBS         = -llapack -lblas -L$(MQCLIB)
F03Flags     = 
RunF         = pgfortran -DPGI -i8 -r8 -Mallocatable=03
#RunF         = pgfortran -i8 -r8
#
#
# The 'all' rule.
#
all: nio.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
%.mod: %.F03
	$(RunF) $(LIBS) $(Prof) -I$(MQCMODS) -c $*.F03

%.o: %.f90
	$(RunF) -I$(MQCMODS) -c $*.f90

%.o: %.F03
	$(RunF) $(F03Flags) -I$(MQCMODS) -c $*.F03

#
# Generic rule for building general executable program (*.exe) from a standard
# f03 source (*.f03) file.
#
%.exe: %.F03 %_mod.F03 %_mod.mod $(MQCLIB)/libmqc.a
	$(RunF) $(LIBS) $(Prof) -I$(MQCMODS) -o $*.exe $*.F03 $*_mod.o $(MQCLIB)/libmqc.a
