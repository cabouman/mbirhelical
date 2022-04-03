#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -q premium
#SBATCH -N 10
#SBATCH -C haswell
#SBATCH -t 16:30:00
#SBATCH -J L096
#SBATCH -A als

module load impi
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so
export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so

export NUM_NODES=10
export NUM_FOCALS=2
export NUM_SOURCES=1
export OMP_NUM_THREADS=32
export OMP_PROC_BIND=true   # new recommendations for hybrid MPI/OpenMP
export OMP_PLACES=threads
export DUAL_ENERGY=0
export DEBUG_MODE=0

name="fulldose_L096_3mm"
forward_model_directory="../data/${name}/forward_model_directory.txt"
info_recon_directory="../data/${name}/info_recon.txt"
prior_directory="../data/${name}/prior_qggmrf.txt"
ce_directory="../data/${name}/ce.txt"
recon_directory="/global/cscratch1/sd/wang1698/L096/fulldose_3mm/recon"

#srun -n $NUM_NODES -c 272 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FORWARD_MODELS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 50
srun -n $NUM_NODES -c 64 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCALS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 70 ${DUAL_ENERGY} ${DEBUG_MODE} ${NUM_SOURCES}







echo
echo " ============================ END DEMO ==========================="

