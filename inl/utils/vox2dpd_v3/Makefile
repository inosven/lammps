CC = g++

####flags
CFLAGS = -O3 -Wall -Wno-unused-result
#CFLAGS = -Og -g -pg -fPIC
CFLAGS += -std=c++0x

####Source files
EXE = vox2dpd_v3
OBJ = vox2dpd_v3.o

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

