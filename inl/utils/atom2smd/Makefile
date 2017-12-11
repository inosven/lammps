CC = g++

####flags
CFLAGS = -std=c++0x -O3
CFLAGS += -fopenmp

####Source files
EXE = atom2smd
OBJ = atom2smd.o

####Executable
VOXCEL2ATOM:$(EXE)
$(EXE):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ)
.C.o:
	$(CC) $(CFLAGS) -c $*.C

#### clean up old builds
clean:
	rm -f *.o
	rm -f $(EXE)
