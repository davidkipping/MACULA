NESTLIBDIR = ..
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
FC=mpif90
OBJFILES = params.o macula.o like.o nestwrap.o main.o

all: macula-nest

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
macula-nest: $(OBJFILES)
	$(FC) -o ../macula-nest $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../macula-nest ../chains/macula-nest-* *~
