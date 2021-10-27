# For Linux or any machines with gcc compiler
CC = cc
CFLAGS = -qopenmp -O3

all: ct fpro clean

clean:
	/bin/rm *.o

OBJ = allocate.o io.o proj.o prepro.o solve.o icd.o randlib.o quicksort.o

ct: ct.o $(OBJ) 
	$(CC) $(CFLAGS) -o ct ct.o $(OBJ) -lm

fpro: fpro.o $(OBJ)
	$(CC) $(CFLAGS) -o fpro fpro.o $(OBJ) -lm

