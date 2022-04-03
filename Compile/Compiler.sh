#module load impi
module load craype-haswell
#module swap craype-haswell craype-mic-knl

cd ../src/
make
