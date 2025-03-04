#!/bin/ksh
set -ex
#
# module load ifort/2019.1.144-GCC-8.2.0-2.31.1
# module load netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
#FC="pgf90"
#FFLAGS="-byteswapio -Kieee -C -Ktrap=fp"
# ON FRAM
#NETCDF_DIR="/cluster/software/netCDF/4.4.1.1-intel-2018a-HDF5-1.8.19/"
NETCDF_DIR="/cluster/software/netCDF-Fortran/4.6.0-iimpi-2022a/"
#
FC="ifort"
FFLAGS=""
LIB="-L$NETCDF_DIR/lib -lnetcdf -lnetcdff" 
LDFLAGS="-g $LIBS"
INC=-I/$NETCDF_DIR/include
#NETCDF_DIR=/cluster/software/netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19/include/ #/work/apps/netcdf/3.6.3-pgi
#HDF5_DIR=/opt/cray/hdf5/1.8.16/pgi/153/
#INC=-I/$(NETCDF_DIR)/include
#LIB="-L$NETCDF_DIR/lib -lnetcdf -lnetcdff -L$HDF5_DIR/lib -lhdf5 -lhdf5_hl"
#LIB="" #"-L$NETCDF_DIR/lib -lnetcdf"

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
