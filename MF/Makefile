CFLAGS=-O3 -fopenmp -Wno-unused-result
PATHS=-I/usr/include -L/usr/lib/x86_64-linux-gnu/
LIBS=-lgsl -lgslcblas -lm

MuFLR: MuFLR.c Makefile AU_fftgrid.h AU_cosmoparam.h AU_cosmofunc.h AU_fastpt_coord.h
	gcc MuFLR.c -o MuFLR $(CFLAGS) $(PATHS) $(LIBS) 

clean:
	$(RM) MuFLR

