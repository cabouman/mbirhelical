#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -N 16
#SBATCH -p batch
#SBATCH -t 08:00:00
#SBATCH -J AAPM_all
#SBATCH -A gen006


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

module load gcc

for i in {130..133}
do

	weight_name="aapm-parameters/dcm_$(printf %03d $i)"
	echo "${weight_name}"

	forward_model_directory="../data/${weight_name}/forward_model_directory.txt"
	info_recon_directory="../data/${weight_name}/info_recon.txt"
	prior_directory="../data/${weight_name}/prior_qggmrf.txt"
	ce_directory="../data/${weight_name}/ce.txt"
	recon_directory="/gpfs/alpine/gen006/scratch/xf9/recon/dcm$(printf %03d $i)/recon"


	srun -N 4 -n 4 -c ${OMP_NUM_THREADS} --cpu_bind=cores  --exclusive ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 100 ${DUAL_ENERGY} ${DEBUG_MODE} ${NUM_SOURCES} &
done

wait

echo
echo " ============================ END DEMO ==========================="

