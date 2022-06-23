module load gcc
# for andes cluster

export TF_DIR=/gpfs/alpine/proj-shared/gen006/muraligm/software/tensorflow_2_8_0_cpp
export TF_INCLUDE_DIR=$TF_DIR/include
export TF_LIBRARY_DIR=$TF_DIR/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TF_LIBRARY_DIR


cd ../src/
make
