#!/bin/bash
#
# Job name:
#SBATCH --job-name=HresGrid
# project
#SBATCH --account=nn2345k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# Max memory usage per core (GB):
#SBATCH --mem-per-cpu=32GB
# Wall clock limit:
#SBATCH --time=00:20:00
#SBATCH --partition=bigmem
#
set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error
#
module --quiet purge
module load ifort/2019.1.144-GCC-8.2.0-2.31.1
module load netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
#
workdir=${USERWORK}/create_Hres_grid
cd $workdir
ln -s depth_regridded.nc depth.nc
./micom_grid
# END
