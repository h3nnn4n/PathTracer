CC=g++
C_FLAGS=$(CFLAGS) -Wall -fopenmp
SOURCES=trig.cpp path.cpp vector.cpp color.cpp magic.cpp plane.cpp sphere.cpp
BIN=magic

all:
	$(CC) $(C_FLAGS) $(SOURCES) -o $(BIN)
