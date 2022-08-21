#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -N 16
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -t 00:30:00
#SBATCH -J AAPM_all

module load impi
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so
export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/usr/common/software/libfabric/1.5.0/gnu/lib/libfabric.so

export NUM_NODES=$SLURM_JOB_NUM_NODES

echo "$NUM_NODES"

export NUM_FOCAL_SPOTS=1
export NUM_SOURCES=1

export OMP_NUM_THREADS=32
#export OMP_NUM_THREADS=68
export OMP_PROC_BIND=true   # new recommendations for hybrid MPI/OpenMP
export OMP_PLACES=threads
export DUAL_ENERGY=0
export DEBUG_MODE=0


#srun -n $NUM_NODES -c 272 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 70 ${DUAL_ENERGY} ${DEBUG_MODE} ${NUM_SOURCES}

for i in {144..147}
do

	weight_name="aapm-parameters/dcm_$(printf %03d $i)"
	echo "${weight_name}"

	forward_model_directory="../data/${weight_name}/forward_model_directory.txt"
	info_recon_directory="../data/${weight_name}/info_recon.txt"
	prior_directory="../data/${weight_name}/prior_qggmrf.txt"
	ce_directory="../data/${weight_name}/ce.txt"
	recon_directory="/global/cscratch1/sd/wang1698/AAPM_2022/recon/dcm$(printf %03d $i)/recon"


	srun -N 4 -n 4 -c 64 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 50 ${DUAL_ENERGY} ${DEBUG_MODE} ${NUM_SOURCES} &
done

wait

echo
echo " ============================ END DEMO ==========================="

