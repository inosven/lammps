CC = g++

####flags
CFLAGS = -O3 -Wall -Wno-unused-result
CFLAGS += -std=c++0x

####Source files
EXE = voxel2dpd
OBJ = voxel2dpd.o

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

