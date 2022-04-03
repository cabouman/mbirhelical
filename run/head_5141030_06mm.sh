#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -q premium
#SBATCH -N 20
#SBATCH -C haswell
#SBATCH -t 16:30:00
#SBATCH -J head
#SBATCH -L SCRATCH

module load impi
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so
export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so

export NUM_NODES=20
export NUM_FOCAL_SPOTS=4
export NUM_SOURCES=1
export OMP_NUM_THREADS=32
export OMP_PROC_BIND=true   # new recommendations for hybrid MPI/OpenMP
export OMP_PLACES=threads
export DUAL_ENERGY=0
export DEBUG_MODE=0



weight_name="head_5141030/head_scan_06mm"
echo "${weight_name}"

forward_model_directory="../data/${weight_name}/forward_model_directory.txt"
info_recon_directory="../data/${weight_name}/info_recon.txt"
prior_directory="../data/${weight_name}/prior_qggmrf.txt"
ce_directory="../data/${weight_name}/ce.txt"
recon_directory="/global/cscratch1/sd/wang1698/head_scan/5141030/MBIR_full_dose_06mm/recon"

#srun -n $NUM_NODES -c 272 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FORWARD_MODELS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 50
srun -n $NUM_NODES -c 64 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 50 ${DUAL_ENERGY} ${DEBUG_MODE} ${NUM_SOURCES}


echo
echo " ============================ END DEMO ==========================="

