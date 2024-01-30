CC=icpx
CFLAGS = -O3 -fiopenmp -fopenmp-targets=spir64 -fopenmp-offload-mandatory -g
#CC=g++
#CFLAGS= -O3 -Wall -Wextra -pedantic -g -lm -fopenmp
#CC=nvc++
#CFLAGS= -O3 -g -lm -mp=gpu
SRCDIR=src
BINDIR=bin
SOURCES=$(wildcard $(SRCDIR)/*.cpp)
OBJECTS=$(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(SOURCES))
EXECUTABLE1=$(BINDIR)/project
EXECUTABLE2=$(BINDIR)/wrt

all: $(EXECUTABLE2) #$(EXECUTABLE1)

$(EXECUTABLE1): $(OBJECTS)
	$(CC) $(CFLAGS) $(EXECUTABLE1).cpp -o $(EXECUTABLE1)_cpp.x $(OBJECTS)

$(EXECUTABLE2): $(OBJECTS)
	$(CC) $(CFLAGS) $(EXECUTABLE2).cpp -o $(EXECUTABLE2)_cpp.x $(OBJECTS)

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE1)_cpp.x $(EXECUTABLE2)_cpp.x
