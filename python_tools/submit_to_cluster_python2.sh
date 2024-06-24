#!/bin/bash
#
# Job name:
#SBATCH --job-name=HresGrid
# project
#SBATCH --account=nn2345k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# Max memory usage per core (GB):
#SBATCH --mem-per-cpu=16GB
# Wall clock limit:
#SBATCH --time=08:00:00
#SBATCH --partition=bigmem
#
set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error
#
#module --quiet purge
#module load ifort/2019.1.144-GCC-8.2.0-2.31.1
#module load netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
#
workdir=${USERWORK}/create_Hres_grid2022_v2
cd $workdir
source /cluster/home/anu074/miniconda3/bin/activate python3env
python remap_inicon.py
# END
