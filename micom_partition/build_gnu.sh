#!/bin/ksh
set -ex
FC="gfortran"
FFLAGS="-fconvert=swap"
LDFLAGS=""
NETCDF_DIR=$CRAY_NETCDF_DIR/gnu/46
HDF5_DIR=$CRAY_HDF5_DIR/gnu/46
INC="-I$NETCDF_DIR/include"
LIB="-L$NETCDF_DIR/lib -lnetcdf -lnetcdff -L$HDF5_DIR/lib -lhdf5 -lhdf5_hl"

$FC $FFLAGS $INC -c types.f90
$FC $FFLAGS $INC -c ncutils.f90
$FC $FFLAGS $INC -c mod_xc.F
$FC $FFLAGS $INC -c mod_za.F
$FC $FFLAGS $INC -c wtime.F
$FC $FFLAGS $INC -c zh.F
$FC $FFLAGS $INC -c partit.f
$FC $FFLAGS $INC -c topo_ppm.f

$FC types.o ncutils.o mod_xc.o wtime.o mod_za.o partit.o zh.o $LDFLAGS $LIB -o partit
$FC types.o ncutils.o mod_xc.o wtime.o mod_za.o topo_ppm.o zh.o $LDFLAGS $LIB -o topo_ppm
