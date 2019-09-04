#
# Makefile for wave_3d_main.c RUN: typing make in this directory
#
#CFLAGS= -O3 -Wall -std=c99
CC=gcc
#CC=icc
# CFLAGS= -O3 -Wall -ggdb3
CFLAGS= -O3 -Wall -m64 -std=c99 #-openmp
# $(CC) $(CFLAGS) 

all: ff_exec ff_main.o adisco_1d.o inidat.o  integ.o rkc.o ff_eq.o  util.o input.o derivs_1d.o pygraph.o 

ff_exec: ff_main.c first_macro_1d.h structs_1d.h ff_main.o inidat.o integ.o rkc.o adisco_1d.o ff_eq.o  util.o input.o derivs_1d.o 
	$(CC) $(CFLAGS) -L. -L/usr/local/lib -L/usr/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -I/usr/include -I/usr/include/hdf5/serial -ffast-math -finline-functions -o ff_exec ff_main.o inidat.o integ.o rkc.o adisco_1d.o ff_eq.o util.o input.o derivs_1d.o pygraph.o -lm -lhdf5 

ff_main.o: first_macro_1d.h structs_1d.h ff_main.c
	$(CC) $(CFLAGS) -c -o ff_main.o ff_main.c

adisco_1d.o: first_macro_1d.h structs_1d.h adisco_1d.c
	$(CC) $(CFLAGS) -c -o adisco_1d.o adisco_1d.c

inidat.o: first_macro_1d.h structs_1d.h inidat.c
	$(CC) $(CFLAGS) -c -o inidat.o inidat.c

integ.o: first_macro_1d.h structs_1d.h integ.c
	$(CC) $(CFLAGS) -c -finline-functions -o integ.o integ.c

rkc.o: first_macro_1d.h structs_1d.h rkc.c
	$(CC) $(CFLAGS) -c -finline-functions -o rkc.o rkc.c
	
ff_eq.o: ff_eq.c first_macro_1d.h structs_1d.h
	$(CC) $(CFLAGS) -c -finline-functions -o ff_eq.o ff_eq.c

input.o: input.c first_macro_1d.h structs_1d.h
	$(CC) $(CFLAGS) -c -o input.o input.c

pygraph.o: pygraph.c pygraph.h
	$(CC) $(CFLAGS) -I/usr/include/hdf5/serial -c -o pygraph.o pygraph.c



#output.o: output.c first_macro_1d.h structs_1d.h gen_1d.h
#	$(CC) $(CFLAGS) -c -o output.o output.c

derivs_1d.o: derivs_1d.c first_macro_1d.h structs_1d.h derivs_1d.h
	$(CC) $(CFLAGS) -c -ffast-math -finline-functions -o derivs_1d.o derivs_1d.c

#binio.o: binio.c
#	$(CC) $(CFLAGS) -c -DLITTLE -o binio.o binio.c

util.o: util.c
	$(CC) $(CFLAGS) -c -o util.o util.c


clean:
	rm -f *.o mhd_exec tranSVFile core *~

limpia: 
	rm -f Dump/run_*
