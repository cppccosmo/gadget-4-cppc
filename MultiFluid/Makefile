CFLAGS=-O3 -fopenmp -Wno-unused-result
PATHS=-I/usr/include -L/usr/lib/x86_64-linux-gnu/
LIBS=-lgsl -lgslcblas -lm

MF: main.c Makefile utils/fftgrid.h utils/cosmoparam.h utils/cosmofunc.h utils/fastpt_coord.h
	gcc main.c -o MF $(CFLAGS) $(PATHS) $(LIBS) 

clean:
	$(RM) MF

