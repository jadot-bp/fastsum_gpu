CC=icpx
CFLAGS= -g -O2 -fsycl -lm #-g -O3 -Wall -Wextra -pedantic -lm
SRCDIR=src
BINDIR=bin
#SOURCES=$(wildcard $(SRCDIR)/*.cpp)
#OBJECTS=$(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(SOURCES))
SOURCES=
OBJECTS=
EXECUTABLE=$(BINDIR)/wrt

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(EXECUTABLE).cpp -o $(EXECUTABLE)_cpp.x $(OBJECTS)

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)_cpp.x
