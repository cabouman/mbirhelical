#!/bin/sh -l
# FILENAME: submit.sh

#SBATCH -N 4
#SBATCH -t 03:00:00
#SBATCH -J AAPM_134
#SBATCH -A gen150
#SBATCH -p batch


module load gcc


export TF_DIR=/gpfs/alpine/proj-shared/gen006/muraligm/software/tensorflow_2_8_0_cpp
export TF_INCLUDE_DIR=$TF_DIR/include
export TF_LIBRARY_DIR=$TF_DIR/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TF_LIBRARY_DIR



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


weight_name="aapm-parameters/dcm_134"
echo "${weight_name}"

forward_model_directory="../data/${weight_name}/forward_model_directory.txt"
info_recon_directory="../data/${weight_name}/info_recon.txt"
prior_directory="../data/${weight_name}/prior_dl.txt"
ce_directory="../data/${weight_name}/ce.txt"
recon_directory="/gpfs/alpine/gen006/proj-shared/xf9/recon/dcm134/recon"

#srun -n $NUM_NODES -c 272 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 70 ${DUAL_ENERGY} ${DEBUG_MODE} ${NUM_SOURCES}

srun -N $NUM_NODES -n $NUM_NODES -c 32 --cpu_bind=cores ../src/ct ${forward_model_directory} ${NUM_FOCAL_SPOTS} ${info_recon_directory} ${prior_directory} ${ce_directory} ${recon_directory} 70 ${DUAL_ENERGY} ${DEBUG_MODE} ${NUM_SOURCES}


echo
echo " ============================ END DEMO ==========================="

