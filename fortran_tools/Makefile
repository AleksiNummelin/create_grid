# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include Makefile.config


.DEFAULT:
	-touch $@
all: micom_grid
micom_fill.o: ./micom_fill.f90 types.o
	$(FC) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	./micom_fill.f90
micom_grid.o: ./micom_grid.f90 types.o mod_xc.o
	$(FC) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	./micom_grid.f90
mod_xc.o: ./mod_xc.F dimensions.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c -I.	./mod_xc.F
ncerr.o: ./ncerr.f90
	$(FC) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	./ncerr.f90
nctools.o: ./nctools.F mod_xc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(SFFLAGS) $(OTHERFLAGS) -c	./nctools.F
sph_length.o: ./sph_length.f90 types.o
	$(FC) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	./sph_length.f90
sph_polyarea.o: ./sph_polyarea.f90 types.o
	$(FC) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	./sph_polyarea.f90
types.o: ./types.f90
	$(FC) $(FFLAGS) $(SF90FLAGS) $(OTHERFLAGS) -c	./types.f90
SRC = ./sph_length.f90 ./micom_grid.f90 ./types.f90 ./ncerr.f90 ./nctools.F ./micom_fill.f90 ./sph_polyarea.f90 ./mod_xc.F dimensions.h
OBJ = sph_length.o micom_grid.o types.o ncerr.o nctools.o micom_fill.o sph_polyarea.o mod_xc.o
clean: neat
	-rm -f .cppdefs $(OBJ) micom_grid
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
micom_grid: $(OBJ) 
	$(LD) $(OBJ) -o micom_grid  $(LDFLAGS)
