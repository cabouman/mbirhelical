# For Linux or any machines with gcc compiler
# For andes
# To show the required mpi -I and -L directories, use: mpicc --showme
CC = gcc-12
#CFLAGS = -Wall -qopenmp  -no-ansi-alias -restrict -xMIC-AVX512 -O3 -no-prec-div
CFLAGS = -g -Wall -fopenmp -Wno-cast-align -Wrestrict -DDEBUG -fno-tree-dce -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.3.sdk/usr/include -I/usr/local/Cellar/open-mpi/4.1.4_2/include -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.3.sdk/usr/lib -L/usr/local/Cellar/open-mpi/4.1.4_2/lib -L/usr/local/opt/libevent/lib -lmpi
#CFLAGS =  -Wall -fopenmp -Wno-cast-align -Wrestrict -O3 -fno-tree-dce
# TODO:  make versions for debug and no debug, other compilers, etc.


all: ct fpro clean

clean:
	/bin/rm *.o

OBJ = allocate.o io.o proj.o prepro.o solve.o icd.o randlib.o quicksort.o

ct: ct.o $(OBJ) 
	$(CC) $(CFLAGS) -o ct ct.o $(OBJ) -lm
	/bin/mv ct ../bin

fpro: fpro.o $(OBJ)
	$(CC) $(CFLAGS) -o fpro fpro.o $(OBJ) -lm
	/bin/mv fpro ../bin

