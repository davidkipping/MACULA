# Compiler
FC = g95

# Compiler flags
FCFLAGS = -O3

# Programs to build
PROGRAMS = macula-plot

# "make" builds all
all: $(PROGRAMS)

# macula build rules
macula-plot.o: macula.o
macula-plot: macula.o

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS) macula-*.dat macula-*.pdf .DS_Store
