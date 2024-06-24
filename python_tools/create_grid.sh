#!/bin/bash
#
# BEFORE STARTING YOU NEED TO
# 
# 1) CREATE A LOCAL PYTHON ENVIRONMENT, for example
#
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# sh Miniconda3-latest-Linux-x86_64.sh
# conda env create -f environment.yml
# 
# you need xarray, numpy, scipy, netcdf4, scikit-learn (skimage), joblib, gsw
# for plotting one needs matplotib
#
# ################################################
# DEFINE PARAMETERS FOR SUBMISSION
# Job name:
#SBATCH --job-name=HresGrid
# project
#SBATCH --account=nn9560k
#SBATCH --ntasks=1
# remap-depth will run parallel on these cpu's - orig 8
#SBATCH --cpus-per-task=8
# Max memory usage per core (GB):
#SBATCH --mem-per-cpu=48GB
# Wall clock limit:
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem
#
# #################################################
set -o errexit  # Exit the script on any error
#
# CREATE WORK DIRECTORY
workdir=${USERWORK}/create_DOTPaleoGrid_34MA_1
mkdir -p $workdir
cd $workdir
# make here a git checkout
#
# copy files
#cp /cluster/shared/noresm/inputdata/ocn/micom/tnx0.125v1/submission_files/* $workdir
#cp /cluster/home/anu074/make_micom_grid_clean/submission_files/* $workdir
#
# CREATE THE LON-LAT INFORMATION
cd python_tools
# this is for old fram
# source /cluster/home/anu074/miniconda3/bin/activate python3env
# this is for nird
conda activate /projects/NS9252K/users/anu074/KeyCLIM
python tripolar.py
conda deactivate
#module load MATLAB/2018b
#/cluster/software/MATLAB/2018b/bin/matlab -nosplash -nodisplay -nojvm -nodesktop -r "run('tripolar.m'); quit;" > matlab_out.txt
#
# CREATE THE FIRST GRID.NC FILE
#
cd ../fortran_tools/
cp ../python_tools/dimensions.h .
# this is for old FRAM
# module --force purge
# module load StdEnv
# module load ifort/2019.1.144-GCC-8.2.0-2.31.1
# module load netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
module purge
module load netCDF-Fortran/4.6.0-iimpi-2022a
#
make clean
make
#
./micom_grid
#
#module unload ifort/2019.1.144-GCC-8.2.0-2.31.1
#module unload netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
module purge
# USE THE GRID FILE INFORMATION TO REMAP THE DEPTH
# ALSO FILL SMALL INLETS AND LAKES
# this is on FRAM
# cd ..
# sbatch submit_to_cluster_python.sh
# sbatch submit_to_cluster_python1.sh
# this on NIRD
cd ../python_tools/
conda activate /projects/NS9252K/users/anu074/KeyCLIM
python remap_depth_xesmf.py
python modify_bathy.py
conda deactivate
#source /cluster/home/anu074/miniconda3/bin/activate python3env
#python remap_depth.py
#python modify_bathy.py
#conda deactivate
#
cd ../grid_out/
ln -s depth_python_xesmf_regridded.nc depth.nc
# RUN THE GRID CREATIION SECOND TIME TO UPDATE DEPTH AND MASKS
#module --force purge
#module load StdEnv
#module load ifort/2019.1.144-GCC-8.2.0-2.31.1
#module load netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
cd ../fortran_tools/
module load netCDF-Fortran/4.6.0-iimpi-2022a
./micom_grid
module unload netCDF-Fortran/4.6.0-iimpi-2022a
#module unload ifort/2019.1.144-GCC-8.2.0-2.31.1
#module unload netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
# CREATE INITIAL CONDITIONS
#source /cluster/home/anu074/miniconda3/bin/activate python3env
cd ../python_tools/
conda activate /projects/NS9252K/users/anu074/KeyCLIM
python remap_inicon.py
conda deactivate
#
module --force purge
module load StdEnv
module load netCDF-Fortran/4.4.4-intel-2018a-HDF5-1.8.19
module load CDO/1.9.3-intel-2018a
#nccopy -k '64-bit offset' /cluster/shared/noresm/inputdata/ocn/micom/tnx0.125v1/inicon_nearest_s2d.nc /cluster/shared/noresm/inputdata/ocn/micom/tnx0.125v1/inicon_200728.nc
#nccopy -k 'NC_64BIT_DATA' /cluster/work/users/anu074/create_Hres_grid2022/inicon_nearest_s2d.nc /cluster/work/users/anu074/create_Hres_grid2022/inicon_tnx0.125_20220913.n
nccopy -k 'cdf5' esmf_grid_tnx0.25_DOTPaleo.nc esmf_grid_tnx0.25_DOTPaleo_cdf5.nc /projects/NS9874K/noresm/grids/inicon_DOTPaleo025/control/inicon_nearest_s2d.nc /projects/NS9874K/noresm/grids/inicon_DOTPaleo025/control/inicon_tnx0.25_DOTPaleo34Ma_control_15102023.nc
#
# END
