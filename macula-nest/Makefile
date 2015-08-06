FC = mpif90
CC = mpicc
CXX = mpiCC
FFLAGS += -cpp -O3 -DMPI -fno-range-check 
CFLAGS += -O3 -DMPI

LAPACKLIB = -llapack

NESTLIBDIR = ./

export FC CC CXX FFLAGS CFLAGS LAPACKLIB

 
AR = ar r  
LINKLIB = ld -shared
 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

default: libnest3.a

all: libnest3.a macula-nest
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 

macula-nest:
	make -C macula-files

clean: 
	-rm $(NESTLIBDIR)/libnest3.*  *.o *.mod
	
cleanall: clean_exec clean clean_macula-nest

clean_exec:
	-rm macula-nest

clean_macula-nest:
	make -C macula-files clean
