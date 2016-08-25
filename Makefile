CC=g++
CFLAGS=-c -Wall --std=c++11 -O0
CFLAGSL=-lpthread -O0
LDFLAGS=

SOURCES=trig.cpp path.cpp vector.cpp color.cpp magic.cpp plane.cpp sphere.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=magic

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(CFLAGSL)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(CFLAGSL)

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(EXECUTABLE)
	$(RM) *.svg
