# For Linux or any machines with gcc compiler
GCC = icc
CFLAGS = -fopenmp -O3

all: ct fpro clean

clean:
	/bin/rm *.o

OBJ = allocate.o io.o proj.o prepro.o solve.o icd.o randlib.o quicksort.o

ct: ct.o $(OBJ) 
	$(GCC) $(CFLAGS) -o ct ct.o $(OBJ) -lm

fpro: fpro.o $(OBJ)
	$(GCC) $(CFLAGS) -o fpro fpro.o $(OBJ) -lm

#%.o:%.c
#	$(GCC) -c $<
