#!/bin/bash
# Begin LSF directives
#BSUB -P gen006
#BSUB -J Grad_test
#BSUB -o Grad_test.o%J
#BSUB -W 2:00
#BSUB -nnodes 5
#BSUB -alloc_flags "nvme smt4"
#BSUB -q debug
# End LSF directives and begin shell commands

nnodes=$(cat ${LSB_DJOB_HOSTFILE} | sort | uniq | grep -v login | grep -v batch | wc -l)

module load open-ce/1.4.0-py39-0
conda activate /gpfs/alpine/gen006/proj-shared/xf9/my_env
# grab nodecount
nodes=($(cat ${LSB_DJOB_HOSTFILE} | sort | uniq | grep -v login | grep -v batch))
nnodes=${#nodes[@]}


source export_DDP_envvars.sh
export TORCH_CUDA_ARCH_LIST=7.0 
export CC=gcc

jsrun --smpiargs="-disable_gpu_hooks" -n $nnodes -a 6 -c 42 -g 6 -r 1\
    --bind=proportional-packed:7 --launch_distribution=packed \
     python -u sinowght_preprocess.py 
