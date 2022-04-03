#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -C haswell
#SBATCH -t 8:30:00
#SBATCH -J L067
#SBATCH -L SCRATCH
#SBATCH -A als

module load impi
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so
export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so

export NUM_FORWARD_MODELS=2
export NUM_NODES=2
export OMP_NUM_THREADS=32
export OMP_PROC_BIND=true   # new recommendations for hybrid MPI/OpenMP
export OMP_PLACES=threads

weight_name="square_root"
echo "${weight_name}"
export DUAL_ENERGY=0
export DEBUG_MODE=0

forward_model_directory="../data/L067_${weight_name}_weight/forward_model_directory.txt"
info_recon_directory="../data/L067_${weight_name}_weight/info_recon.txt"
prior_directory="../data/L067_${weight_name}_weight/prior_qggmrf.txt"
ce_directory="../data/L067_${weight_name}_weight/ce.txt"
recon_directory="/global/cscratch1/sd/wang1698/L067/${weight_name}_weight/recon"

#srun -n $NUM_NODES -c 272 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FORWARD_MODELS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 50
srun -n $NUM_NODES -c 64 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FORWARD_MODELS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 40 ${DUAL_ENERGY} ${DEBUG_MODE}


echo
echo " ============================ END DEMO ==========================="

