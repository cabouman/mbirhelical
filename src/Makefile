# For Linux or any machines with gcc compiler
CC = cc
#CFLAGS = -Wall -qopenmp  -no-ansi-alias -restrict -xMIC-AVX512 -O3 -no-prec-div
CFLAGS =  -Wall -qopenmp -no-ansi-alias -restrict -xCORE-AVX2 -O3 -no-prec-div


all: ct clean

clean:
	/bin/rm *.o

OBJ = allocate.o io.o proj.o prepro.o solve.o icd.o randlib.o quicksort.o

ct: ct.o $(OBJ) 
	$(CC) $(CFLAGS) -o ct ct.o $(OBJ) -lm
