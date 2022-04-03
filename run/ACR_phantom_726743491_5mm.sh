#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -q regular
#SBATCH -N 4
#SBATCH -C haswell
#SBATCH -t 08:30:00
#SBATCH -J ACR_head_5
#SBATCH -L SCRATCH
#SBATCH -A als

module load impi
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so
export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so

export NUM_NODES=4
export NUM_FOCAL_SPOTS=4
export OMP_NUM_THREADS=32
#export OMP_NUM_THREADS=68
export OMP_PROC_BIND=true   # new recommendations for hybrid MPI/OpenMP
export OMP_PLACES=threads

weight_name="ACR_phantom_726743491_5mm"
echo "${weight_name}"

forward_model_directory="../data/${weight_name}/forward_model_directory.txt"
info_recon_directory="../data/${weight_name}/info_recon.txt"
prior_directory="../data/${weight_name}/prior_qggmrf.txt"
ce_directory="../data/${weight_name}/ce.txt"
recon_directory="/global/cscratch1/sd/wang1698/ACR_phantoms/726743491/5mm/recon"

#srun -n $NUM_NODES -c 272 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 70
srun -n $NUM_NODES -c 64 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 70


echo
echo " ============================ END DEMO ==========================="
