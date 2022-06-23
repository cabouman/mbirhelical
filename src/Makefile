# For Linux or any machines with gcc compiler
# For andes
CC = mpicc
#CFLAGS = -Wall -qopenmp  -no-ansi-alias -restrict -xMIC-AVX512 -O3 -no-prec-div
CFLAGS =  -Wall -fopenmp -Wno-cast-align -Wrestrict -O3 -fno-tree-dce -I${TF_INCLUDE_DIR} -L${TF_LIBRARY_DIR} -ltensorflow 



all: ct clean

clean:
	/bin/rm *.o

OBJ = allocate.o io.o proj.o prepro.o solve.o icd.o randlib.o quicksort.o

ct: ct.o $(OBJ) 
	$(CC) $(CFLAGS) -o ct ct.o $(OBJ) -lm
