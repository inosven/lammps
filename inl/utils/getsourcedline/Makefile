CC = g++

####flags
CFLAGS = -std=c++0x -O3

####Source files
EXE = getsourcedline
OBJ = getsourcedline.o

####Executable
GETSOURCEDLINE:$(EXE)
$(EXE):$(OBJ)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ)
.C.o:
	$(CC) $(CFLAGS) -c $*.C

#### clean up old builds
clean:
	rm -f *.o
	rm -f $(EXE)
