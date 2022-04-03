#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -q premium
#SBATCH -N 4
#SBATCH -C haswell
#SBATCH -t 04:30:00
#SBATCH -J lung_5
#SBATCH -L SCRATCH
#SBATCH -A als

module load impi
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so
export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so

export NUM_NODES=4
export NUM_FOCAL_SPOTS=2
export OMP_NUM_THREADS=32
export OMP_PROC_BIND=true   # new recommendations for hybrid MPI/OpenMP
export OMP_PLACES=threads
export DUAL_ENERGY=0
export DEBUG_MODE=0


weight_name="lung_cancer"
echo "${weight_name}"

forward_model_directory="../data/${weight_name}_5mm/forward_model_directory.txt"
info_recon_directory="../data/${weight_name}_5mm/info_recon.txt"
prior_directory="../data/${weight_name}_5mm/prior_qggmrf.txt"
ce_directory="../data/${weight_name}_5mm/ce.txt"
recon_directory="/global/cscratch1/sd/wang1698/lung_cancer/MBIR_full_dose_5mm/recon"

#srun -n $NUM_NODES -c 272 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FORWARD_MODELS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 50
srun -n $NUM_NODES -c 64 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 120 ${DUAL_ENERGY} ${DEBUG_MODE}



echo
echo " ============================ END DEMO ==========================="

