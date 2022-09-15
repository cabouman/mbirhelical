# For Linux or any machines with gcc compiler
# To show the required mpi -I and -L directories, use: mpicc --showme
CC = /usr/local/bin/gcc-12
#CFLAGS = -Wall -qopenmp  -no-ansi-alias -restrict -xMIC-AVX512 -O3 -no-prec-div
CFLAGS =  -Wall -fopenmp -Wno-cast-align -Wrestrict -D DEBUG -fno-tree-dce -I/usr/local/Cellar/open-mpi/4.1.4_2/include -L/usr/local/Cellar/open-mpi/4.1.4_2/lib -L/usr/local/opt/libevent/lib -lmpi
#CFLAGS =  -Wall -fopenmp -Wno-cast-align -Wrestrict -O3 -fno-tree-dce

all: ct fpro clean

clean:
	/bin/rm *.o

OBJ = allocate.o io.o proj.o prepro.o solve.o icd.o randlib.o quicksort.o

ct: ct.o $(OBJ)
	$(CC) $(CFLAGS) -o ct ct.o $(OBJ) -lm

fpro: fpro.o $(OBJ)
	$(CC) $(CFLAGS) -o ct ct.o $(OBJ) -lm


