#!/bin/bash
#
# Job name:
#SBATCH --job-name=PaleoGrid
# project
#SBATCH --account=nn9869k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# Max memory usage per core (GB):
#SBATCH --mem-per-cpu=8GB
# Wall clock limit:
#SBATCH --time=08:00:00
#SBATCH --partition=bigmem
#
set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error
#
workdir=${USERWORK}/create_DOTPaleoGrid_34MA_1/python_tools/
cd $workdir
source /cluster/home/anu074/miniconda3/bin/activate python3env
python remap_depth_xesmf.py
# END
